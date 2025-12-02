module FDatasets

export DSForcing, DSResults, last_date

using Dates

using NCDatasets
using Oceananigans
using ClimaOcean.DataWrangling: Metadatum, metadata_path

import Oceananigans: location
import ClimaOcean.DataWrangling:
    metadata_filename,
    default_download_directory,
    all_dates,
    first_date,
    last_date,
    dataset_variable_name,
    download_dataset,
    longitude_interfaces,
    latitude_interfaces,
    z_interfaces,
    reversed_vertical_axis,
    inpainted_metadata_path

struct DSForcing
    metadata_filename::String
    default_download_directory::String
    reversed_vertical_axis::Any
    longitude_interfaces::Any
    latitude_interfaces::Any
    size::Any
    all_dates::Any
    first_date::Any
    last_date::Any
    z_interfaces::Any
end

function DSForcing(metadata_filename, default_download_directory)
    reversed_vertical_axis = false
    filepath = joinpath(default_download_directory, metadata_filename)
    ds = NCDataset(filepath)
    longitude_interfaces = (ds["Nx"][1], ds["Nx"][end])
    latitude_interfaces = (ds["Ny"][1], ds["Ny"][end])
    size = (ds.dim["Nx"], ds.dim["Ny"], ds.dim["Nz"])
    all_dates = ds["time"][:]
    first_date = ds["time"][1]
    last_date = ds["time"][end]
    z_interfaces = ds["Nz_faces"][:]
    close(ds)

    return DSForcing(
        metadata_filename,
        default_download_directory,
        reversed_vertical_axis,
        longitude_interfaces,
        latitude_interfaces,
        size,
        all_dates,
        first_date,
        last_date,
        z_interfaces,
    )
end  # function

struct DSResults
    metadata_filename::String
    default_download_directory::String
    reversed_vertical_axis::Any
    longitude_interfaces::Any
    latitude_interfaces::Any
    size::Any
    all_dates::Any
    first_date::Any
    last_date::Any
    z_interfaces::Any
end

function DSResults(metadata_filename, default_download_directory; date_time)
    reversed_vertical_axis = false
    filepath = joinpath(default_download_directory, metadata_filename)
    ds = NCDataset(filepath)
    longitude_interfaces = (ds["λ_faa"][1], ds["λ_faa"][end])
    latitude_interfaces = (ds["φ_afa"][1], ds["φ_afa"][end])
    array_size = Dict(:tracers => size(ds["T"])[1:3], :u_velocity => size(ds["u"])[1:3], :v_velocity => size(ds["v"])[1:3])
    all_dates = date_time .+ Second.(ds["time"][:])
    first_date = all_dates[1]
    last_date = all_dates[end]
    z_interfaces = ds["z_aaf"][:]
    close(ds)

    return DSResults(
        metadata_filename,
        default_download_directory,
        reversed_vertical_axis,
        longitude_interfaces,
        latitude_interfaces,
        array_size,
        all_dates,
        first_date,
        last_date,
        z_interfaces,
    )
end  # function

Variable_names = Dict(:temperature => "T", :salinity => "S", :u_velocity => "u", :v_velocity => "v")
Variable_location = Dict(
    :temperature => (Center, Center, Center),
    :salinity => (Center, Center, Center),
    :free_surface => (Center, Center, Nothing),
    :sea_ice_thickness => (Center, Center, Nothing),
    :sea_ice_concentration => (Center, Center, Nothing),
    :net_heat_flux => (Center, Center, Nothing),
    :u_velocity => (Face, Center, Center),
    :v_velocity => (Center, Face, Center),
    :sensible_heat_flux => (Center, Center, Nothing),
    :latent_heat_flux => (Center, Center, Nothing),
    :net_longwave => (Center, Center, Nothing),
    :downwelling_shortwave => (Center, Center, Nothing),
    :downwelling_longwave => (Center, Center, Nothing),
)

const MetadatumForcing = Metadatum{<:DSForcing,<:Any,<:Any}
const MetadatumResults = Metadatum{<:DSResults,<:Any,<:Any}

metadata_filename(metadata::Union{MetadatumForcing, MetadatumResults}) = metadata.dataset.metadata_filename
default_download_directory(ds::Union{DSResults, DSForcing}) = ds.default_download_directory
reversed_vertical_axis(ds::Union{DSResults, DSForcing}) = ds.reversed_vertical_axis
longitude_interfaces(ds::Union{DSResults, DSForcing}) = ds.longitude_interfaces
latitude_interfaces(ds::Union{DSResults, DSForcing}) = ds.latitude_interfaces
Base.size(ds::DSForcing) = ds.size
Base.size(ds::DSForcing, variable) = size(ds)
function Base.size(ds::DSResults, variable) 
    if variable == :u_velocity
        return ds.size[:u_velocity]
    elseif variable == :v_velocity
        return ds.size[:v_velocity]
    else
        return ds.size[:tracers]
    end
end

all_dates(ds::Union{DSResults, DSForcing}, args...) = ds.all_dates
first_date(ds::Union{DSResults, DSForcing}, args...) = ds.first_date
last_date(ds::Union{DSResults, DSForcing}, args...) = ds.last_date

z_interfaces(metadata::Union{MetadatumResults, MetadatumForcing}) = metadata.dataset.z_interfaces

function dataset_variable_name(metadata::Union{MetadatumResults, MetadatumForcing}) 
    if haskey(Variable_names, metadata.name) 
        return Variable_names[metadata.name]
    end
    return metadata.name
end

function location(metadata::Union{MetadatumResults, MetadatumForcing})
    if haskey(Variable_location, metadata.name) 
        return Variable_location[metadata.name]
    end
    # assume a tracer
    @info string("Assuming ", metadata.name, " is at (Center, Center, Center).")
    return (Center, Center, Center)
end

function download_dataset(metadata::Union{MetadatumResults, MetadatumForcing})
    filepath = metadata_path(metadata)
    return filepath
end

function inpainted_metadata_filename(metadata::Union{MetadatumResults, MetadatumForcing})
    original_filename = metadata_filename(metadata)
    without_extension = original_filename[1:end-3]
    return without_extension * "_inpainted.jld2"
end

inpainted_metadata_path(metadata::Union{MetadatumResults, MetadatumForcing}) = joinpath(metadata.dir, inpainted_metadata_filename(metadata))

end  # module
