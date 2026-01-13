module NORA3

export NORA3PrescribedAtmosphere

using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: λnodes, φnodes, on_architecture
using Oceananigans.Fields: interpolate!
using Oceananigans.OutputReaders: Cyclical, TotallyInMemory, AbstractInMemoryBackend, FlavorOfFTS, time_indices
using ClimaOcean
using ClimaOcean.OceanSeaIceModels: PrescribedAtmosphere, TwoBandDownwellingRadiation
using Adapt
using NCDatasets
using JLD2
using Dates

import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: new_backend, update_field_time_series!

struct MultiYearNORA3
    metadata_filename::String
    default_download_directory::String
    size::Any
    all_dates::Any
end

function MultiYearNORA3(metadata_filename::String, default_download_directory::String)
    filepath = joinpath(default_download_directory, metadata_filename)
    ds = NCDataset(filepath)
    array_size = size(ds["x_wind_10m"])[1:2]
    all_dates = ds["time"][:]
    close(ds)

    return MultiYearNORA3(metadata_filename, default_download_directory, array_size, all_dates)
end  # function

const NORA3Metadata{D} = Metadata{<:MultiYearNORA3,D}

Base.size(metadata::NORA3Metadata) = metadata.dataset.size

is_three_dimensional(data::NORA3Metadata) = false
location(::JRA55Metadata) = (Center, Center, Center)

all_dates(ds::MultiYearNORA3, name) = ds.all_dates
all_dates(ds::MultiYearNORA3) = ds.all_dates
first_date(ds::MultiYearNORA3) = first(all_dates(ds))
last_date(ds::MultiYearNORA3) = last(all_dates(ds))

# this is field data series stuff to get data during simulation
struct NORA3NetCDFBackend{M} <: AbstractInMemoryBackend{Int}
    start::Int
    length::Int
    metadata::M
end

Adapt.adapt_structure(to, b::NORA3NetCDFBackend) = NORA3NetCDFBackend(b.start, b.length, nothing)

NORA3NetCDFBackend(length, metadata::Metadata) = NORA3NetCDFBackend(1, length, metadata)
NORA3NetCDFBackend(start::Integer, length::Integer) = NORA3NetCDFBackend(start, length, nothing)
NORA3NetCDFBackend(length) = NORA3NetCDFBackend(1, length, nothing)

Base.length(backend::NORA3NetCDFBackend) = backend.length
Base.summary(backend::NORA3NetCDFBackend) = string("NORA3NetCDFBackend(", backend.start, ", ", backend.length, ")")

const NORA3NetCDFFTS = FlavorOfFTS{<:Any,<:Any,<:Any,<:Any,<:NORA3NetCDFBackend}

new_backend(b::NORA3NetCDFBackend, start, length) = NORA3NetCDFBackend(start, length, b.metadata)

function NORA3FieldTimeSeries(variable_name::Symbol, architecture, FT; dataset, start_date, end_date, kw...)

    native_dates = all_dates(dataset, variable_name)
    dates = compute_native_date_range(native_dates, start_date, end_date)
    metadata = Metadata(variable_name; dataset, dates, dataset.default_download_directory)

    return NORA3FieldTimeSeries(metadata, architecture, FT; kw...)
end

function NORA3FieldTimeSeries(
    metadata::NORA3Metadata, architecture, FT;
    latitude = nothing,
    longitude = nothing,
    backend = InMemory(),
    time_indexing = Cyclical(),
)

    # Cannot use `TotallyInMemory` backend with MultiYearJRA55 dataset
    if metadata.dataset isa MultiYearJRA55 && backend isa TotallyInMemory
        msg = string("The `InMemory` backend is not supported for the MultiYearJRA55 dataset.")
        throw(ArgumentError(msg))
    end

    # First thing: we download the dataset!
    download_dataset(metadata)

    # Regularize the backend in case of `JRA55NetCDFBackend`
    if backend isa JRA55NetCDFBackend
        if backend.metadata isa Nothing
            backend = JRA55NetCDFBackend(backend.length, metadata)
        end

        if backend.length > length(metadata)
            backend = JRA55NetCDFBackend(backend.start, length(metadata), metadata)
        end
    end

    # Unpack metadata details
    dataset = metadata.dataset
    name = metadata.name
    time_indices = JRA55_time_indices(dataset, metadata.dates, name)

    # Change the metadata to reflect the actual time indices
    dates = all_dates(dataset, name)[time_indices]
    metadata = Metadata(metadata.name; dataset = metadata.dataset, dates, dir = metadata.dir)

    shortname = dataset_variable_name(metadata)
    variable_name = metadata.name

    filepath = metadata_path(metadata) # Might be multiple paths!!!
    filepath = filepath isa AbstractArray ? first(filepath) : filepath

    # OnDisk backends do not support time interpolation!
    # Disallow OnDisk for JRA55 dataset loading
    if ((backend isa InMemory) && !isnothing(backend.length)) || backend isa OnDisk
        msg = string(
            "We cannot load the JRA55 dataset with a $(backend) backend. Use `InMemory()` or `JRA55NetCDFBackend(N)` instead.",
        )
        throw(ArgumentError(msg))
    end

    if !(variable_name ∈ JRA55_variable_names)
        variable_strs = Tuple("  - :$name \n" for name in JRA55_variable_names)
        variables_msg = prod(variable_strs)

        msg = string(
            "The variable :$variable_name is not provided by the JRA55-do dataset!",
            '\n',
            "The variables provided by the JRA55-do dataset are:",
            '\n',
            variables_msg,
        )

        throw(ArgumentError(msg))
    end

    # Record some important user decisions
    totally_in_memory = backend isa TotallyInMemory

    # Determine default time indices
    if totally_in_memory
        # In this case, the whole time series is in memory.
        # Either the time series is short, or we are doing a limited-area
        # simulation, like in a single column. So, we conservatively
        # set a default `time_indices = 1:2`.
        time_indices_in_memory = time_indices
        native_fts_architecture = architecture
    else
        # In this case, part or all of the time series will be stored in a file.
        # Note: if the user has provided a grid, we will have to preprocess the
        # .nc JRA55 data into a .jld2 file. In this case, `time_indices` refers
        # to the time_indices that we will preprocess;
        # by default we choose all of them. The architecture is only the
        # architecture used for preprocessing, which typically will be CPU()
        # even if we would like the final FieldTimeSeries on the GPU.
        time_indices_in_memory = 1:length(backend)
        native_fts_architecture = architecture
    end

    ds = Dataset(filepath)

    # Note that each file should have the variables
    #   - ds["time"]:     time coordinate
    #   - ds["lon"]:      longitude at the location of the variable
    #   - ds["lat"]:      latitude at the location of the variable
    #   - ds["lon_bnds"]: bounding longitudes between which variables are averaged
    #   - ds["lat_bnds"]: bounding latitudes between which variables are averaged
    #   - ds[shortname]: the variable data

    # Nodes at the variable location
    λc = ds["lon"][:]
    φc = ds["lat"][:]

    # Interfaces for the "native" JRA55 grid
    λn = Array(ds["lon_bnds"][1, :])
    φn = Array(ds["lat_bnds"][1, :])

    # The netCDF coordinates lon_bnds and lat_bnds do not include
    # the last interfaces, so we push them here.
    push!(φn, 90)
    push!(λn, λn[1] + 360)

    i₁, i₂, j₁, j₂, TX = compute_bounding_indices(longitude, latitude, nothing, Center, Center, λc, φc)

    λr = λn[i₁:i₂+1]
    φr = φn[j₁:j₂+1]
    Nrx = length(λr) - 1
    Nry = length(φr) - 1
    close(ds)

    N = (Nrx, Nry)
    H = min.(N, (3, 3))

    JRA55_native_grid = LatitudeLongitudeGrid(
        native_fts_architecture,
        FT;
        halo = H,
        size = N,
        longitude = λr,
        latitude = φr,
        topology = (TX, Bounded, Flat),
    )

    boundary_conditions = FieldBoundaryConditions(JRA55_native_grid, (Center(), Center(), nothing))
    start_time = first_date(metadata.dataset, metadata.name)
    times = native_times(metadata; start_time)

    if backend isa JRA55NetCDFBackend
        fts = FieldTimeSeries{Center,Center,Nothing}(
            JRA55_native_grid,
            times;
            backend,
            time_indexing,
            boundary_conditions,
            path = filepath,
            name = shortname,
        )

        set!(fts)
        return fts
    else
        fts = FieldTimeSeries{Center,Center,Nothing}(
            JRA55_native_grid,
            times;
            time_indexing,
            backend,
            boundary_conditions,
        )

        # Fill the data in a GPU-friendly manner
        ds = Dataset(filepath)
        data = ds[shortname][i₁:i₂, j₁:j₂, time_indices_in_memory]
        close(ds)

        copyto!(interior(fts, :, :, 1, :), data)
        fill_halo_regions!(fts)

        return fts
    end
end

function NORA3PrescribedAtmosphere(
    architecture = CPU(),
    FT = Float32;
    dataset = MultiYearNORA3("NORA3.nc", joinpath(homedir(), "FjordSim_data", "NORA3")),
    start_date = first_date(dataset),
    end_date = last_date(dataset),
    backend = NORA3NetCDFBackend(10),
    time_indexing = Cyclical(),
    surface_layer_height = 10,  # meters
    other_kw...,
)

    kw = (; time_indexing, backend, start_date, end_date, dataset)
    kw = merge(kw, other_kw)

    ua = NORA3FieldTimeSeries(:eastward_velocity, architecture, FT; kw...)
    va = NORA3FieldTimeSeries(:northward_velocity, architecture, FT; kw...)
    Ta = NORA3FieldTimeSeries(:temperature, architecture, FT; kw...)
    qa = NORA3FieldTimeSeries(:specific_humidity, architecture, FT; kw...)
    pa = NORA3FieldTimeSeries(:sea_level_pressure, architecture, FT; kw...)
    Fra = NORA3FieldTimeSeries(:rain_freshwater_flux, architecture, FT; kw...)
    Fsn = NORA3FieldTimeSeries(:snow_freshwater_flux, architecture, FT; kw...)
    Ql = NORA3FieldTimeSeries(:downwelling_longwave_radiation, architecture, FT; kw...)
    Qs = NORA3FieldTimeSeries(:downwelling_shortwave_radiation, architecture, FT; kw...)

    freshwater_flux = (rain = Fra, snow = Fsn)

    times = ua.times
    grid = ua.grid

    velocities = (u = ua, v = va)

    tracers = (T = Ta, q = qa)

    pressure = pa

    downwelling_radiation = TwoBandDownwellingRadiation(shortwave = Qs, longwave = Ql)

    FT = eltype(ua)
    surface_layer_height = convert(FT, surface_layer_height)

    auxiliary_freshwater_flux = nothing
    atmosphere = PrescribedAtmosphere(
        grid,
        times;
        velocities,
        freshwater_flux,
        auxiliary_freshwater_flux,
        tracers,
        downwelling_radiation,
        surface_layer_height,
        pressure,
    )

    return atmosphere
end # function

end # module