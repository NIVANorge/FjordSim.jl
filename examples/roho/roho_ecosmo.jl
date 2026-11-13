# `roho.jl` with the NIVA ECOSMO biogeochemistry (NivaEcosmo.jl): its 17 tracers beside T and S over
# 2023, on the same grid, physics, rivers and NorKyst-800m boundaries.
#
# # Running it
#
# This file needs NivaEcosmo, which FjordSim's own environment does not have, so it runs in the
# environment beside it, `examples/roho/Project.toml`:
#
#   julia --project=examples/roho -e 'using Pkg; Pkg.instantiate()'
#   julia --project=examples/roho -m FjordSim <step> --config examples/roho/roho_ecosmo.jl
#
# The NorKyst downloads, the bathymetry and the atmosphere are `roho.jl`'s, read from its data
# directory, so run its steps first. The initial state, the boundaries and the rivers carry the
# ECOSMO tracers, so this file prepares its own into `roho_ecosmo`: `download_initial_state`,
# `prepare_initial_state`, `add_rivers`, `download_boundary_data` and `prepare_boundary_data`.
#
# On the GCP VMs it runs in an image built with that environment (`gcp/README.md`, "A setup with
# its own environment"), with `roho`'s data staged beside its own. The environment also carries YAML,
# so one image runs NivaEcosmo's `fjordsim_model/moro.jl` too.
#
#   gcp/fjordsim-gcp run --config examples/roho/roho_ecosmo.jl --stage roho --steps ...
#
# # Where the tracer values come from
#
# **Sea** (open boundaries and initial state): TOPAZ5-ECOSMO (Copernicus Marine
# ARCTIC_ANALYSISFORECAST_BGC_002_004, MET Norway; HYCOM with ECOSMO II at 6.25 km, daily means),
# whose water off Rogaland reaches 266 m. It gives nitrate, phosphate, silicate, oxygen and DIC
# directly, and phytoplankton carbon, chlorophyll and zooplankton carbon as totals this file splits
# between ECOSMO's groups. Ammonium, detritus, labile DOM, opal, calcite and alkalinity it lacks, and
# they are seasonal profiles from the literature below. TOPAZ resolves no fjord, so the initial
# fjord water takes the nearest coastal water TOPAZ has.
#
# **Rivers**: every river carries the same seasonal concentrations, from the 2018 monthly chemistry
# of Suldalslågen, Vikedalselva, Bjerkreimselva, Frafjordelva, Espedalselva and Jørpelandsåna [A].
#
# **Every value marked "assumed" has no measurement behind it.** Replace them with Vannmiljø or the
# river monitoring programme's 2023 data (Elveovervåkingsprogrammet, RID), or with ICES casts on the
# Rogaland shelf, when those are at hand.
#
# [A] Hellen and Hindar (2019), Kalking i laksevassdrag skadet av sur nedbør. Tiltaksovervåking i
#     2018, Miljødirektoratet M-1566, appendix B.
# [B] ØKOKYST delprogram Nordsjøen sør, årsrapport 2019, Miljødirektoratet M-1604: nutrients 0-30 m
#     and phytoplankton groups in Hidlefjorden, Hjelmelandsfjorden and Jøsenfjorden, 2016-2019.
# [C] Monitoring ocean acidification in Norwegian seas in 2019, Miljødirektoratet M-1735.
# [D] Frigstad et al. (2020), Influence of riverine input on Norwegian coastal systems,
#     Front. Mar. Sci. 7:332.
# [E] Alvarez-Fernandez and Riegman (2014), J. Sea Res. 91:35-44: North Sea C:Chl of 12-119 g/g.

using FjordSim
using NivaEcosmo
using NivaEcosmo: ECOSMO_SINKING_TRACERS
# Loads NivaEcosmo's FjordSim extension, which writes the tracers into the prepared files.
using NCDatasets: NCDatasets
using Oceananigans: WENO
using Oceananigans.Units: day, hour
using Dates: DateTime, dayofyear

base = include("roho.jl")

# --- Units ------------------------------------------------------------------------------------

# ECOSMO carries nutrients, plankton and organic matter in mgC m⁻³ through its Redfield ratios.
const MGC_PER_MMOL_N = NivaEcosmo.REDFIELD_C_N * NivaEcosmo.CARBON_MOLAR_MASS
const MGC_PER_MMOL_P = NivaEcosmo.REDFIELD_C_P * NivaEcosmo.CARBON_MOLAR_MASS
const MGC_PER_MMOL_SI = NivaEcosmo.REDFIELD_C_SI * NivaEcosmo.CARBON_MOLAR_MASS
const MGC_PER_MMOL_C = NivaEcosmo.CARBON_MOLAR_MASS
# Seawater density, kg L⁻¹, from µmol kg⁻¹ to mmol m⁻³.
const SEAWATER_DENSITY = 1.026

"""
    seasonal(anchors)

`date -> value`, linear in day of year between the `(day, value)` `anchors` and periodic over the year.
"""
function seasonal(anchors)
    days = [first(anchors[end]) - 365; first.(anchors); first(anchors[1]) + 365]
    values = [last(anchors[end]); last.(anchors); last(anchors[1])]
    return function (date)
        day = dayofyear(date)
        upper = findfirst(>(day), days)
        weight = (day - days[upper-1]) / (days[upper] - days[upper-1])
        return (1 - weight) * values[upper-1] + weight * values[upper]
    end
end

# `surface` at the surface, decaying to `deep` with e-folding depth `scale` (m); `z` negative down.
surface_to_deep(z, surface, deep, scale) = deep + (surface - deep) * exp(z / scale)

# --- The sea: TOPAZ5-ECOSMO and the profiles it lacks ------------------------------------------

# Diatoms' share of phytoplankton carbon: diatoms (mostly Skeletonema) dominate the spring bloom,
# February-May, and flagellates the rest of the year, in the Ryfylke fjords [B]. The shares are
# assumed, since [B] counts cells rather than carbon.
const DIATOM_FRACTION = seasonal([(25, 0.25), (40, 0.75), (140, 0.75), (155, 0.25)])

# Carbon to chlorophyll, g g⁻¹, splitting TOPAZ's chlorophyll between the groups. Within the North
# Sea's 12-119 [E], diatoms lower than flagellates; the two values are assumed.
const DIATOM_CARBON_TO_CHLOROPHYLL = 40.0
const FLAGELLATE_CARBON_TO_CHLOROPHYLL = 80.0

# Microzooplankton's share of zooplankton carbon: half in winter, a third through the growing
# season. Assumed: no regional measurement was found.
const MICROZOOPLANKTON_FRACTION = seasonal([(80, 0.5), (95, 0.35), (150, 0.35), (165, 0.3), (320, 0.3), (340, 0.5)])

# Ammonium, µM: 1.1-1.8 at 0-10 m in the Ryfylke fjords and nearly flat through the year [B]; the
# deep value and its 50 m scale are assumed.
const AMMONIUM_SURFACE = seasonal([(105, 1.3), (150, 1.2), (275, 1.2), (320, 1.3)])
const AMMONIUM_DEEP = 0.2

# Dissolved organic carbon, µM: 2.5-3 mg L⁻¹ (210-250 µM) on the Norwegian coast from spring to
# autumn [D], taken lower at the open boundary; the winter and deep values and the 100 m scale are
# assumed. ECOSMO's `dom` is the labile pool only, an assumed 5 % of it.
const ORGANIC_CARBON_SURFACE = seasonal([(90, 120.0), (150, 200.0), (320, 200.0), (350, 120.0)])
const ORGANIC_CARBON_DEEP = 50.0
const LABILE_FRACTION = 0.05

# Detritus, mgC m⁻³: low in winter and summer, peaks in spring and autumn, as in the central North
# Sea; the values and the 150 m scale are assumed.
const DETRITUS_SURFACE = seasonal([(30, 20.0), (100, 60.0), (150, 30.0), (240, 30.0), (290, 60.0), (340, 20.0)])
const DETRITUS_DEEP = 5.0

# Opal, µM Si, following the diatom bloom; calcite, mgC m⁻³, small. Both assumed, with 50 m scales.
const OPAL_SURFACE = seasonal([(25, 0.2), (40, 1.0), (140, 1.0), (155, 0.2)])
const OPAL_DEEP = 0.1
const CALCITE_SURFACE = 0.5
const CALCITE_DEEP = 0.1

"""
    alkalinity(S)

Total alkalinity, mmol m⁻³, as the mixing line between the rivers' 70 µmol kg⁻¹ [A] and Atlantic
water's 2315 µmol kg⁻¹ at S = 35 [C], capped there.
"""
alkalinity(S) = SEAWATER_DENSITY * min(70 + 64.1 * S, 2320)

"""
    inorganic_carbon(dic, S)

TOPAZ's DIC, mmol m⁻³, diluted towards the rivers' 150 below the salinity of the coastal water it
represents, 34. TOPAZ resolves no fjord, so a fresh cell near a river mouth would otherwise carry
seawater's DIC beside river water's alkalinity.
"""
inorganic_carbon(dic, S) = 150 + (dic - 150) * clamp(S / 34, 0, 1)

"""
    chlorophyll_share(fla, dia)

Diatoms' share of the chlorophyll of `fla` and `dia` mgC m⁻³ at their carbon-to-chlorophyll ratios.
"""
function chlorophyll_share(fla, dia)
    flagellates = fla / FLAGELLATE_CARBON_TO_CHLOROPHYLL
    diatoms = dia / DIATOM_CARBON_TO_CHLOROPHYLL
    return diatoms / max(flagellates + diatoms, eps())
end

phytoplankton(source) = MGC_PER_MMOL_C * max(source.phyc, 0)
zooplankton(source) = MGC_PER_MMOL_C * max(source.zooc, 0)
flagellates(date, source) = (1 - DIATOM_FRACTION(date)) * phytoplankton(source)
diatoms(date, source) = DIATOM_FRACTION(date) * phytoplankton(source)
diatom_chlorophyll(date, source) =
    max(source.chl, 0) * chlorophyll_share(flagellates(date, source), diatoms(date, source))

# Every ECOSMO tracer at the open boundaries and in the initial state, from TOPAZ (`source.no3`, ...,
# its units) and from the base file's salinity (`source.S`).
const SEA = (
    no3 = (z, date, source) -> MGC_PER_MMOL_N * max(source.no3, 0),
    nh4 = (z, date, source) -> MGC_PER_MMOL_N * surface_to_deep(z, AMMONIUM_SURFACE(date), AMMONIUM_DEEP, 50),
    pho = (z, date, source) -> MGC_PER_MMOL_P * max(source.po4, 0),
    sil = (z, date, source) -> MGC_PER_MMOL_SI * max(source.si, 0),
    oxy = (z, date, source) -> max(source.o2, 0),
    fla = (z, date, source) -> flagellates(date, source),
    dia = (z, date, source) -> diatoms(date, source),
    flachl = (z, date, source) -> max(source.chl, 0) - diatom_chlorophyll(date, source),
    diachl = (z, date, source) -> diatom_chlorophyll(date, source),
    microzoo = (z, date, source) -> MICROZOOPLANKTON_FRACTION(date) * zooplankton(source),
    mesozoo = (z, date, source) -> (1 - MICROZOOPLANKTON_FRACTION(date)) * zooplankton(source),
    det = (z, date, source) -> surface_to_deep(z, DETRITUS_SURFACE(date), DETRITUS_DEEP, 150),
    opa = (z, date, source) -> MGC_PER_MMOL_SI * surface_to_deep(z, OPAL_SURFACE(date), OPAL_DEEP, 50),
    cal = (z, date, source) -> surface_to_deep(z, CALCITE_SURFACE, CALCITE_DEEP, 50),
    dom = (z, date, source) -> LABILE_FRACTION * MGC_PER_MMOL_C *
                               surface_to_deep(z, ORGANIC_CARBON_SURFACE(date), ORGANIC_CARBON_DEEP, 100),
    dic = (z, date, source) -> inorganic_carbon(1000 * source.dissic, source.S),
    alk = (z, date, source) -> alkalinity(source.S),
)

const TOPAZ_PARAMETERS = ["no3", "po4", "si", "o2", "dissic", "chl", "phyc", "zooc"]

# --- The rivers ---------------------------------------------------------------------------------

# 2018 monthly samples [A]: nitrate 74-346 µg N L⁻¹ across the stations, highest in late winter;
# silicate 0.7-1.5 mg SiO₂ L⁻¹, highest in winter; alkalinity 0.04-0.12 mmol L⁻¹; TOC 2.35 mg L⁻¹,
# taken as DOC. Phosphorus sits at the 2 µg L⁻¹ detection limit, so phosphate is assumed at a third
# of it. Ammonium, DIC (from the alkalinity at pH 6.3), detritus and oxygen at saturation are assumed.
# The plankton, opal and calcite enter at zero.
cosine(mean, amplitude, peak_day) = date -> mean + amplitude * cos(2π * (dayofyear(date) - peak_day) / 365)

const RIVERS = (
    no3 = date -> MGC_PER_MMOL_N * cosine(12.5, 2.5, 70)(date),
    nh4 = date -> MGC_PER_MMOL_N * 0.35,
    pho = date -> MGC_PER_MMOL_P * 0.03,
    sil = date -> MGC_PER_MMOL_SI * cosine(17.5, 6.0, 40)(date),
    oxy = cosine(357.5, 42.5, 30),
    det = date -> 300.0,
    dom = date -> LABILE_FRACTION * MGC_PER_MMOL_C * 190.0,
    dic = date -> 150.0,
    alk = date -> 70.0,
)

# --- Light, wind and CO₂ ------------------------------------------------------------------------

# Days since midnight on 1 January, as ECOSMO's `start_day_of_year` counts them.
const START_DAY_OF_YEAR = dayofyear(start_date) - 1
const SOLAR_CONSTANT = 1361.0
# Assumed, as on `moro.jl`: NORA3's shortwave would be better, but is not reachable from here.
const CLEARNESS_INDEX = 0.5
# 10 m wind for the gas exchange, m s⁻¹, assumed; CO₂ in air, ppm, the 2023 global mean.
const WIND_SPEED = 7.0
const ATMOSPHERIC_PCO2 = 420.0

"""
    surface_shortwave(λ, φ, t)

Daily-mean downwelling shortwave below the surface, W m⁻², at `t` seconds after `start_date`: the
top-of-atmosphere daily mean (Cooper's declination) times `CLEARNESS_INDEX`. Runs in a kernel.
"""
@inline function surface_shortwave(λ, φ, t)
    day_of_year = START_DAY_OF_YEAR + t / 86400
    declination = deg2rad(23.44) * sin(2π * (day_of_year + 284) / 365)
    latitude = deg2rad(φ)
    hour_angle = acos(clamp(-tan(latitude) * tan(declination), -1, 1))
    eccentricity = 1 + 0.033 * cos(2π * day_of_year / 365)
    insolation = SOLAR_CONSTANT / π * eccentricity *
                 (hour_angle * sin(latitude) * sin(declination) + cos(latitude) * cos(declination) * sin(hour_angle))
    return CLEARNESS_INDEX * max(0, insolation)
end

# --- The setup ----------------------------------------------------------------------------------

ecosmo_configuration = fabm_configuration(joinpath(pkgdir(NivaEcosmo), "data", "fabm.yaml"), FT)

# This file's own data directory; `data_root` is `roho.jl`'s, whose downloads it reads.
ecosmo_root = fjord_data_root("roho_ecosmo")

# A copy of one of `roho.jl`'s components with `fields` changed.
function changed(config; fields...)
    copy = deepcopy(config)
    for (field, value) in pairs(fields)
        setproperty!(copy, field, value)
    end
    return copy
end

physics = base.model_config

FjordConfig(
    grid_config = base.grid_config,
    bathymetry_config = base.bathymetry_config,
    initial_state_config = ecosmo_initial_state(
        changed(base.initial_state_config; data_root = ecosmo_root,
                output_directory = joinpath(data_root, base.initial_state_config.output_directory));
        sources = TopazInitialStateConfig(
            data_root = ecosmo_root,
            output_directory = "topaz_initial",
            parameters = TOPAZ_PARAMETERS,
            start_date = start_date,
        ),
        tracers = SEA,
    ),
    river_config = ecosmo_rivers(
        changed(base.river_config; data_root = ecosmo_root),
        RIVERS,
    ),
    boundary_data_config = ecosmo_boundaries(
        changed(base.boundary_data_config; data_root = ecosmo_root,
                output_directory = joinpath(data_root, base.boundary_data_config.output_directory));
        sources = TopazBoundaryDataConfig(
            data_root = ecosmo_root,
            output_directory = "topaz_daily",
            open_edges = base.boundary_data_config.open_edges,
            parameters = TOPAZ_PARAMETERS,
            start_date = start_date,
            stop_date = stop_date,
        ),
        tracers = SEA,
    ),
    atmosphere_config = base.atmosphere_config,
    model_config = CoupledHydrostaticModel(
        buoyancy = physics.buoyancy,
        closure = physics.closure,
        # The ECOSMO tracers on plain WENO, ECOSMO's state repair clipping their undershoots, but
        # bounded for the three that sink: on `moro.jl` plain WENO drained the cell above the seabed
        # below zero, the repair clipped it back, and the mass created piled up in the bottom cell.
        tracer_advection = merge(
            physics.tracer_advection,
            NamedTuple{ECOSMO_TRACERS}(map(_ -> WENO(FT), ECOSMO_TRACERS)),
            NamedTuple{ECOSMO_SINKING_TRACERS}(map(_ -> WENO(FT; bounds = FT.((0, 1e6))), ECOSMO_SINKING_TRACERS)),
        ),
        momentum_advection = physics.momentum_advection,
        # The whole list: the river sources, the open-boundary conditions and the initial conditions
        # are decided from it before the model exists.
        tracers = (:T, :S, ECOSMO_TRACERS...),
        coriolis = physics.coriolis,
        sea_ice = physics.sea_ice,
        radiative_forcing = physics.radiative_forcing,
        free_surface = physics.free_surface,
        extra_kwargs = physics.extra_kwargs,
    ),
    biogeochemistry_config = ecosmo_biogeochemistry_config(
        configuration = ecosmo_configuration,
        surface_shortwave = surface_shortwave,
        # `roho.jl`'s bottom drag, which the sediment's resuspension stress should share.
        bottom_drag_coefficient = 0.003,
        start_day_of_year = START_DAY_OF_YEAR,
    ),
    boundary_condition_configs = (
        base.boundary_condition_configs...,
        ecosmo_gas_exchange(
            wind_speed = WIND_SPEED,
            atmospheric_pCO2 = ATMOSPHERIC_PCO2,
            pH_scale = ecosmo_configuration.pH_scale,
        ),
    ),
    run_config = RunConfig(
        results_root = fjord_results_root("roho_ecosmo"),
        architecture = base.run_config.architecture,
        writers = (
            # `roho.jl`'s writers but its surface one, which the ECOSMO surface writer replaces.
            filter(writer -> !(writer isa SurfaceSnapshotWriter), base.run_config.writers)...,
            # Daily, a new file every month. One record of the 17 tracers is about 220 MB, so the
            # year is about 80 GB before compression.
            SnapshotWriter(
                name = :ecosmo,
                output_file = "snapshots_ecosmo.nc",
                variables = ECOSMO_TRACERS,
                interval = 1day,
                split_times = monthly_splits,
                overwrite_existing = true,
            ),
            # Hourly, a new file every month: the surface of T, S, u, v, η and every ECOSMO field.
            ecosmo_surface_writer(
                name = :surface,
                output_file = "snapshots_surface.jld2",
                interval = 1hour,
                split_times = monthly_splits,
                overwrite_existing = true,
            ),
        ),
        callbacks = base.run_config.callbacks,
        time_stepping = base.run_config.time_stepping,
        # NorKyst's T, S, u and v with the ECOSMO tracers appended (`initial_state_config`).
        initial_conditions = FromInitialState(),
        start_date = start_date,
        stop_date = stop_date,
        loops = 1,
        pickup = false,
    ),
)
