# The `examples/oslofjorden.jl` setup with the simplest OceanBioME.jl biogeochemistry bolted on: a
# four-compartment nitrogen NPZD model (nutrient, phytoplankton, zooplankton, detritus) in
# mmol N m⁻³.
#
# Nothing in `src/` changes. The point of this file is to show the three ways an out-of-tree config
# extends FjordSim, in increasing order of reach:
#
#   1. **Reuse.** `oslofjorden.jl` is included, and the pieces that need no change — the grid, the
#      bathymetry, the atmosphere — are passed straight through.
#   2. **A new config subtype with a new hook method.** `NPZDRivers` and `NPZDBoundaries` subtype an
#      `Abstract*Config`, hold the built-in config they extend, and delegate every hook but the one
#      or two they mean to change; `NPZDBiogeochemistry` implements `AbstractBiogeochemistryConfig`'s
#      one hook. That is the documented way to add a source (`docs/adding-a-source.md`) and it is all
#      rivers and the model need.
#   3. **A new method on a *pipeline* generic.** `prepare_boundary_data` is a generic function
#      dispatching on the config type, so `NPZDBoundaries` adds a method that runs the built-in
#      pipeline and then *appends* the biogeochemical variables to the file it wrote.
#
# (2) and (3) are what buy everything else for free. Once `boundaries_npzd.nc` carries `south_N` and
# the river series carry `N`, two pieces of FjordSim that were written to be name-driven pick the new
# tracers up with no further code:
#
#   * `river_forcing` builds a source for every model tracer the river file carries, so river nitrate
#     enters with the same water, through the same columns, as river temperature and salinity. P, Z
#     and D, which no river names, enter at the value of the water the river flows into.
#   * `open_tracer_boundary_conditions` radiates *every* model tracer towards its own prepared
#     series, so `OpenLateralBoundaryFromData` needs no change at all — it would in fact **error**
#     without `south_N`, which is why the boundary append is not optional.
#
# N, P, Z and D start from the same profiles, evaluated at `start_date` and given to
# `initial_conditions` as functions of depth.
#
# The prepared files are new names under the shared `data_root`, so the downloads are reused and
# `oslofjorden.jl`'s own `rivers_nve.nc` and `boundaries.nc` are left untouched:
#
#   julia --project -m FjordSim prepare_bathymetry  --config examples/oslofjorden.jl   # shared, if not run
#   julia --project -m FjordSim download_boundary_data --config examples/oslofjorden.jl   # shared
#   julia --project -m FjordSim download_atmosphere --config examples/oslofjorden.jl   # shared
#   julia --project -m FjordSim prepare_atmosphere  --config examples/oslofjorden.jl   # shared
#   julia --project -m FjordSim add_rivers          --config examples/oslofjorden_npzd.jl
#   julia --project -m FjordSim prepare_boundary_data  --config examples/oslofjorden_npzd.jl
#   julia --project -m FjordSim run_simulation      --config examples/oslofjorden_npzd.jl
#
# `add_rivers` needs `NVE_API_KEY`, exactly as `oslofjorden.jl` does.

using FjordSim
using FjordSim: BoundaryData, Rivers, Simulations
using FjordSim.Configs:
    AbstractBiogeochemistryConfig, AbstractBoundaryDataConfig, AbstractRiverConfig
using Oceananigans
using Oceananigans.Units
using OceanBioME
using NCDatasets
using Dates: dayofyear

# --- The literature ---------------------------------------------------------------------------
#
# Every number below is from NIVA report L.NR. 7626-2021, *Overvåking av Ytre Oslofjord 2019-2023.
# Tilførsler og undersøkelser i vannmassene i 2020* — the monitoring year this run simulates, so
# the profiles and the forcing describe the same twelve months. Where a value is an assumption
# rather than a measurement it says so.

# Station OF-1 (outer Oslofjord, 440 m) on 10 March 2020: the pre-bloom profile, the one visit of
# the year when nitrate is closest to conservative, so it is the shape to hang a seasonal cycle on.
# Reported as NO₃+NO₂ in µg N/L; divided by 14.007 to reach mmol N m⁻³.
#
# The surface node is not OF-1's — the report samples that station from 50 m down — but station
# LA-1 at 2 m, which read 94 and 95 µg N/L on 15 January and 12 February 2020.
const NITRATE_PROFILE_DEPTHS = [0.0, 50.0, 75.0, 100.0, 125.0, 150.0, 200.0, 250.0, 300.0, 400.0, 440.0]
const NITRATE_PROFILE_WINTER = [6.7, 6.4, 7.9, 8.6, 9.3, 8.6, 8.6, 8.6, 9.3, 12.8, 13.6]

# Summer surface nitrate, from the same LA-1 visits in June, July and August 2020, every one of
# which came back below the 1 µg N/L detection limit. Not zero, because a nutrient a biogeochemical
# model divides by should not be.
const NITRATE_SUMMER_SURFACE = 0.05

# The drawdown is a growing-season *window* rather than a peak, because that is what the data shows:
# LA-1 read 94 µg N/L at 2 m in January, then below detection on every visit from 18 June to
# 15 September, then 95 µg N/L again on 10 November. A Gaussian centred on midsummer cannot hold a
# four-month floor and still come back up by November — fitted to both ends it sat 15x above the
# observed August surface — so the shape is a product of two logistics: drawdown switching on
# through the spring bloom, and off again when autumn mixing re-supplies the surface.
#
# The onset is the Oslofjord spring bloom, which runs from late February into April; the end is
# mid-October. Both widths are the fortnight-scale of those transitions.
const BLOOM_ONSET_YEARDAY = 90.0
const BLOOM_ONSET_WIDTH = 18.0
const BLOOM_END_YEARDAY = 290.0
const BLOOM_END_WIDTH = 12.0

# How far down the drawdown reaches, as an e-folding scale. 50 m is what the OF-1 series asks for
# and it is deeper than a photic-zone argument would give: on 18 August 2020 the station read
# 49 µg N/L at 50 m against 115 at 100 m, and on 16 September 1.2 µg N/L at 50 m — so the summer
# deficit is still a third of the winter value two euphotic depths down.
const NITRATE_DEPLETION_SCALE = 50.0

# Summer chlorophyll a at 2 m, June-August 2020 (report Table 7): 1.0-4.7 µg/L across the sixteen
# stations, 2.5 µg/L being representative of the outer ones this domain is judged on. Converted at
# a carbon-to-chlorophyll ratio of 50 g C / g Chl and Redfield C:N = 106:16 by mole (5.68 g C / g N),
# so 2.5 µg Chl/L is 22 µg N/L, or 1.6 mmol N m⁻³.
const PHYTOPLANKTON_SUMMER = 1.6

# The report measures chlorophyll only from June to September, so the winter minimum is an
# assumption, not a measurement: 3 % of the summer value, which is the usual overwintering
# inoculum a spin-up needs to have something to grow from.
const PHYTOPLANKTON_WINTER = 0.05

# Both assumptions. No Oslofjord mesozooplankton or detrital nitrogen series was found for 2020, and
# these are the conventional NPZD spin-up fractions — a grazer stock well below its prey, and a
# detritus pool small enough that the first weeks of the run set it rather than this file.
const ZOOPLANKTON_FRACTION = 0.3
const DETRITUS_FRACTION = 0.2

# The euphotic e-folding scale the plankton profiles use. Shallower than the nitrate one because
# biomass tracks light, not the mixed layer.
const EUPHOTIC_SCALE = 25.0

# How well the nitrate profile above actually does, scored against all 23 OF-1 and LA-1 nitrate
# measurements the report publishes for 2020: mean bias +0.9, RMSE 2.2 mmol N m⁻³ on a field that
# spans 0-13. The two ends that matter are close — 6.65 against 6.71 at the surface in January,
# which is the instant the initial conditions are evaluated at, and 8.4 / 8.6 / 12.8 against 8.6 / 8.6 / 12.9
# at 100, 200 and 400 m in the March profile.
#
# Every large residual is subsurface and late-season: 4.1 against 0.1 at 50 m in September, 7.8
# against 2.2 at 100 m the same day. Those are not a seasonal cycle this shape is missing — OF-1's
# 50 m nitrate runs 90, 6.6, 11, 49, 1.2, 22 µg N/L through the year, which is Skagerrak water
# masses arriving and leaving, not a smooth drawdown. Reproducing them is the model's job, driven by
# the open boundary; an analytic profile that tried would be fitting advection with a curve.

# Flow-weighted mean total nitrogen, in mmol N m⁻³, for the rivers the report gauges (its Table 4,
# 2019): the published annual load divided by the published annual water volume, times the nitrate
# fraction below. Glomma is 16 007 t N on 68 724 × 10³ m³ d⁻¹, which is 638 µg TotN/L.
#
# Numedalslågen and Skienselva are outside this domain — both reach the sea west of 10.2°E — and are
# here only because they set `RIVER_NITRATE_DEFAULT`.
const RIVER_TOTAL_NITROGEN = Dict(
    "Glomma (Osterelva)" => 638.0,
    "Glomma (Vesterelva)" => 638.0,
    "Drammenselva" => 463.0,
)

# What fraction of a river's total nitrogen is dissolved inorganic, which is what NPZD's `N` is.
# Measured in the fjord rather than in the rivers, at the two Drammensfjord stations in 2020
# (report Table 6): NO₃+NO₂ over TOTN was 157/293 and 177/347 in summer, and 270/410 at D-2 on
# 14 January. 0.6 is the middle of that, and it is the weakest number on this page.
const RIVER_NITRATE_FRACTION = 0.6

# Every other river in the domain — the nineteen `NVERiversConfig` discovers and this file never
# names — gets the four gauged rivers' combined flow-weighted mean, 506 µg TotN/L: 24 849 t N on
# 134 516 × 10³ m³ d⁻¹. A domain-mean concentration is a much better guess for an ungauged stream
# than Glomma's, and it keeps the load scaling with each stream's own discharge, which is what the
# river forcing already does with volume.
const RIVER_TOTAL_NITROGEN_DEFAULT = 506.0

# µg N/L to mmol N m⁻³, at the atomic weight of nitrogen.
const NITROGEN_MOLAR_MASS = 14.007

# Surface irradiance. `PAR_FRACTION` is the standard photosynthetically active share of shortwave
# energy; `CLEARNESS_INDEX` is southern Norway's annual-mean ratio of surface to top-of-atmosphere
# insolation, from Oslo's ~950 kWh m⁻² yr⁻¹ global horizontal irradiation (108 W m⁻²) against an
# annual-mean top-of-atmosphere ~210 W m⁻² at this latitude.
#
# This is the one place the example fabricates a field FjordSim already has: `prescribed_radiation`
# reads NORA3's downwelling shortwave for the air-sea fluxes, and feeding *that* to the light model
# would be strictly better than a clear-sky formula with a constant cloud factor. It is not done
# here because it would mean reaching inside NumericalEarth's radiation object, which is a bigger
# change than this file is for.
const SOLAR_CONSTANT = 1361.0
const PAR_FRACTION = 0.43
const CLEARNESS_INDEX = 0.5

"""
    seasonal_drawdown(date)

How much of the winter-to-summer nitrate deficit is in place on `date`, from 0 in midwinter to 1
through the productive season.

The one seasonal shape this file has, shared by the nutrient and the plankton profiles so that
biomass appears exactly where nitrate disappears. Sharing it is what makes the plankton profiles a
summer plateau rather than a spring peak, which is as much as the report supports — it measures
chlorophyll only from June to September, and those are the values `PHYTOPLANKTON_SUMMER` is set by.
"""
function seasonal_drawdown(date)
    day = dayofyear(date)
    onset = 1 / (1 + exp(-(day - BLOOM_ONSET_YEARDAY) / BLOOM_ONSET_WIDTH))
    ending = 1 / (1 + exp((day - BLOOM_END_YEARDAY) / BLOOM_END_WIDTH))
    return onset * ending
end

"""
    winter_nitrate(z)

The March 2020 OF-1 nitrate profile at depth `z` (negative down, metres), linearly interpolated
between the reported levels and held constant below the deepest one.
"""
function winter_nitrate(z)
    depth = -z
    depth <= first(NITRATE_PROFILE_DEPTHS) && return first(NITRATE_PROFILE_WINTER)
    depth >= last(NITRATE_PROFILE_DEPTHS) && return last(NITRATE_PROFILE_WINTER)

    upper = searchsortedlast(NITRATE_PROFILE_DEPTHS, depth)
    span = NITRATE_PROFILE_DEPTHS[upper+1] - NITRATE_PROFILE_DEPTHS[upper]
    weight = (depth - NITRATE_PROFILE_DEPTHS[upper]) / span

    return (1 - weight) * NITRATE_PROFILE_WINTER[upper] + weight * NITRATE_PROFILE_WINTER[upper+1]
end

"""
    nitrate_profile(z, date)

Dissolved inorganic nitrogen at depth `z` on `date`, in mmol N m⁻³: the winter profile less a
surface-intensified seasonal deficit.

The deficit is scaled by the *surface* winter value rather than by the local one, so the same
absolute drawdown is removed at every depth the euphotic scale reaches — which is what OF-1 shows.
Floored at the summer surface value so the deep column can never be driven negative by a shallow
argument.
"""
function nitrate_profile(z, date)
    deficit = (first(NITRATE_PROFILE_WINTER) - NITRATE_SUMMER_SURFACE) * seasonal_drawdown(date)
    return max(NITRATE_SUMMER_SURFACE, winter_nitrate(z) - deficit * exp(z / NITRATE_DEPLETION_SCALE))
end

"""
    phytoplankton_profile(z, date)

Phytoplankton nitrogen at depth `z` on `date`, in mmol N m⁻³: the winter stock plus a summer bloom
confined to the euphotic layer, on `seasonal_drawdown`'s own calendar.
"""
phytoplankton_profile(z, date) =
    PHYTOPLANKTON_WINTER +
    (PHYTOPLANKTON_SUMMER - PHYTOPLANKTON_WINTER) * seasonal_drawdown(date) * exp(z / EUPHOTIC_SCALE)

zooplankton_profile(z, date) = ZOOPLANKTON_FRACTION * phytoplankton_profile(z, date)
detritus_profile(z, date) = DETRITUS_FRACTION * phytoplankton_profile(z, date)

# The four profiles, keyed by the tracer name they are written under. One `NamedTuple` because both
# the initial conditions and the boundary append want exactly the same set, evaluated the same way.
const NPZD_PROFILES = (
    N = nitrate_profile,
    P = phytoplankton_profile,
    Z = zooplankton_profile,
    D = detritus_profile,
)

"""
    river_nitrate(name)

Dissolved inorganic nitrogen in the river whose outlet is called `name`, in mmol N m⁻³.

Discovered rivers the report does not gauge fall back to the domain-mean concentration, so an
ungauged stream still carries nitrogen in proportion to its own discharge — the volume, and hence
the load, comes from the river inflow itself.
"""
river_nitrate(name) =
    get(RIVER_TOTAL_NITROGEN, name, RIVER_TOTAL_NITROGEN_DEFAULT) * RIVER_NITRATE_FRACTION /
    NITROGEN_MOLAR_MASS

"""
    surface_par(λ, φ, t)

Daily-mean photosynthetically active radiation at the sea surface, W m⁻², at longitude `λ` and
latitude `φ` (degrees) and `t` seconds into the run.

Clear-sky top-of-atmosphere insolation — Cooper's (1969) declination series and the standard
daily-mean integral over the hour angle — scaled by a constant clearness index and the PAR fraction.
At 59.5°N it gives 110 W m⁻² at the solstice and 6 W m⁻² in January, either side of OceanBioME's
own 100 W m⁻² default.

This is a boundary condition on OceanBioME's `PAR` field, so it runs inside a kernel: `ifelse` and
`clamp` rather than branches, and no allocation. `t` is seconds from the run's `start_date`, which
this file sets to 1 January — so day-of-year is `t / 86400` with no offset, and moving `start_date`
means adding one here.
"""
@inline function surface_par(λ, φ, t)
    day = t / 86400
    declination = deg2rad(23.44) * sin(2π * (day + 284) / 365)
    latitude = deg2rad(φ)

    # Sunrise hour angle. Clamped because |tan φ · tan δ| exceeds one inside the polar circles, and
    # a kernel cannot afford `acos` to leave its domain even where this grid never goes.
    hour_angle = acos(clamp(-tan(latitude) * tan(declination), -1, 1))
    eccentricity = 1 + 0.033 * cos(2π * day / 365)

    insolation =
        SOLAR_CONSTANT / π * eccentricity * (
            hour_angle * sin(latitude) * sin(declination) +
            cos(latitude) * cos(declination) * sin(hour_angle)
        )

    return PAR_FRACTION * CLEARNESS_INDEX * max(0, insolation)
end

# --- Writing the profiles into the prepared files -----------------------------------------------

# Matches FjordSim's own writer settings, so the appended variables compress like the ones
# `prepare_boundary_data` wrote beside them.
const NPZD_DEFLATE_LEVEL = 5

"""
    append_boundary_variables!(filepath, edges, profiles)

Add one `<edge>_<name>` variable per entry of `profiles` and per edge to the prepared boundary file
at `filepath`, on the layout `prepare_boundary_data` writes: the along-edge axis, then depth, then time.

Uniform along the edge and finite everywhere, including under land. The boundary reader fills
non-finite cells from their nearest finite neighbour anyway (`fill_boundary_gaps!`), so a mask here
would be work with no effect — and the exterior state is a property of the water arriving, not of
which model cells happen to be wet.

`boundary_variable_name` rather than string interpolation, because the `<edge>_<name>` convention is
FjordSim's to state and `boundary_series` reads the file back through the same function.
"""
function append_boundary_variables!(filepath, edges, profiles)
    NCDataset(filepath, "a") do ds
        z = ds["Nz"][:]
        dates = ds["time"][:]
        depths = Float32[profile(depth, date) for profile in values(profiles), depth in z, date in dates]

        for edge in edges
            along = edge in (:south, :north) ? "Nx" : "Ny"
            n_along = ds.dim[along]

            for (row, name) in enumerate(keys(profiles))
                variable_name = BoundaryData.boundary_variable_name(edge, String(name))
                @info "Writing $variable_name into $(basename(filepath))"

                variable = defVar(
                    ds, variable_name, Float32, (along, "Nz", "time");
                    deflatelevel = NPZD_DEFLATE_LEVEL,
                    attrib = ["_FillValue" => NaN32],
                )
                variable[:, :, :] = repeat(
                    reshape(depths[row, :, :], 1, length(z), length(dates)), n_along, 1, 1,
                )
            end
        end
    end

    return filepath
end

# --- NPZDBiogeochemistry: a new AbstractBiogeochemistryConfig -----------------------------------

"""
    NPZDBiogeochemistry(; surface_PAR, scale_negatives)

OceanBioME's NPZD, as the knobs `model_biogeochemistry` builds it from.

A config rather than `NPZD(grid; ...)` itself, because `NPZD` needs the grid and this file is
evaluated before there is one — `build_simulation` reads the grid out of the processed bathymetry
and only then calls `model_biogeochemistry`, handing the result to `coupled_model`.
"""
struct NPZDBiogeochemistry{P} <: AbstractBiogeochemistryConfig
    surface_PAR::P
    scale_negatives::Bool
end

NPZDBiogeochemistry(; surface_PAR, scale_negatives) =
    NPZDBiogeochemistry(surface_PAR, Bool(scale_negatives))

"""
    negative_tracer_scaler(scale_negatives)

The `modifiers` an `NPZDBiogeochemistry` asks for: negative scaling over the four biogeochemical
tracers, or `nothing`.

Spelled out rather than left to `NPZD`'s `scale_negatives` keyword, which would include temperature
in the group — see `model_biogeochemistry` below.

`invalid_fill_value = 0` rather than `ScaleNegativeTracers`' own `NaN` default: a cell whose N+P+Z+D
total goes negative is rare but not impossible even away from `T`, and the thinnest partial-bottom
cells this bathymetry has are numerically fragile enough to hit it within the first few hundred
seconds of a run. `NaN` there is not a local error — every neighbouring cell's tracers are `NaN` by
the next advection step, and the whole active domain follows within one output interval. `0` gives up
exact mass conservation at that one cell in exchange for not losing the run.
"""
negative_tracer_scaler(scale_negatives::Bool) =
    scale_negatives ? ScaleNegativeTracers(Tuple(keys(NPZD_PROFILES)); invalid_fill_value = 0) : nothing

"""
    model_biogeochemistry(config::NPZDBiogeochemistry, grid)

OceanBioME's `NPZD` on `grid`. It returns an already-wrapped `Biogeochemistry`; passing it to
`HydrostaticFreeSurfaceModel` as it comes is correct, and wrapping it again is not.

`scale_negatives` is worth having on for a run like this one: WENO on a sharp detritus or
phytoplankton front will undershoot, and a negative nutrient is how these models blow up rather than
merely drift. It cannot be `NPZD`'s own `scale_negatives = true`, though, and that is not a matter of
taste. That keyword builds the modifier from `conserved_tracers(bgc)`, which groups the tracers by
element through `required_biogeochemical_tracers` — and `PhytoZoo`'s list is `(:P, :Z, :T)`, `:T`
being *temperature*, which it requires for its Q₁₀ and which is not a nitrogen pool at all. The
scaler would therefore rescale temperature along with the nutrients on any cell holding a negative
tracer, and set it to `NaN` wherever `N + P + Z + T + D` went negative. A `NaN` temperature is a
`NaN` buoyancy, hence a `NaN` velocity, hence a `NaN` `Δt` — which is exactly how this run used to
die, in `calculate_substeps`, before the first hour was out.

So the modifier is built by `negative_tracer_scaler` instead, over the four tracers that really do
carry the nitrogen, and `scale_negatives` stays at its `false` default. `NPZD_PROFILES` names them, so
the set cannot drift apart from the one the initial conditions and the boundary append use.
"""
Simulations.model_biogeochemistry(config::NPZDBiogeochemistry, grid) = NPZD(
    grid;
    surface_PAR = config.surface_PAR,
    modifiers = negative_tracer_scaler(config.scale_negatives),
)

# --- NPZDRivers: a new AbstractRiverConfig ------------------------------------------------------

"""
    NPZDRivers(base; output_file, plot_file)

`base`'s rivers, carrying nitrate as well as temperature and salinity.

Every hook but `river_series` delegates, so the outlets are still the ones `NVERiversConfig`
discovers from NVE's ELVIS network, still entering through the same columns with the same
discharge. Only the values change, and only by one variable.

The alternative was `NVERiversConfig`'s `constants` field, which already writes one fixed value per
variable at every outlet. It is rejected here for what it cannot say: nitrogen in a Norwegian river
is a catchment property — Glomma drains the agricultural Østlandet at 638 µg TotN/L, Drammenselva
463 — and one number for all twenty-one outlets would throw that away.
"""
struct NPZDRivers{R} <: AbstractRiverConfig
    base::R
    data_root::String
    output_file::String
    plot_file::String
    search_radius::Int
    maximum_inflow_velocity::Float64
end

NPZDRivers(base; output_file, plot_file) = NPZDRivers(
    base,
    base.data_root,
    String(output_file),
    String(plot_file),
    base.search_radius,
    base.maximum_inflow_velocity,
)

Rivers.river_locations(rivers::NPZDRivers) = Rivers.river_locations(rivers.base)
Rivers.river_discharge(rivers::NPZDRivers, times) = Rivers.river_discharge(rivers.base, times)
Rivers.river_search_radius(rivers::NPZDRivers) = Rivers.river_search_radius(rivers.base)
Rivers.river_maximum_inflow_velocity(rivers::NPZDRivers) =
    Rivers.river_maximum_inflow_velocity(rivers.base)
Rivers.download_rivers(target_grid, rivers::NPZDRivers) =
    Rivers.download_rivers(target_grid, rivers.base)

"""
    river_series(rivers::NPZDRivers, times)

`base`'s series plus an `N` row per outlet, in `river_locations` order.

Constant in time, which is a deliberate limit rather than an oversight: NIVA publishes an annual
*load* and an annual water volume, and their ratio is a flow-weighted mean concentration with no
seasonal information in it. A seasonal cycle here would be invention, and the discharge these
concentrations multiply is already the observed daily HydAPI series — so the river nitrogen *flux*
follows the floods even though its concentration does not.

Keyed by outlet *name*, which is the only identifier a `RiverLocation` carries — `vassdragsnr`, what
`oslofjorden.jl` keys its own overrides by, does not survive into one. A name is a weaker key: the
gauged rivers are named only because that setup's overrides state a `name`, and renaming one there
would drop it to `RIVER_TOTAL_NITROGEN_DEFAULT` with nothing to notice. Hence the log line, which
says which outlets were matched and which took the default.
"""
function Rivers.river_series(rivers::NPZDRivers, times)
    series = Rivers.river_series(rivers.base, times)
    locations = Rivers.river_locations(rivers)

    series["N"] = Float32[
        river_nitrate(location.name) for location in locations, _ in times
    ]

    gauged = [location.name for location in locations if haskey(RIVER_TOTAL_NITROGEN, location.name)]
    @info "River nitrate: $(length(gauged)) gauged outlet(s) — $(join(gauged, ", ")) — and " *
          "$(length(locations) - length(gauged)) at the domain mean " *
          "$(round(river_nitrate(""); digits = 1)) mmol N m⁻³"

    return series
end

# --- NPZDBoundaries: a new AbstractBoundaryDataConfig -------------------------------------------

"""
    NPZDBoundaries(base)

`base`'s open-boundary data with the four NPZD tracers appended to the file it prepares.

Not optional. `open_tracer_boundary_conditions` puts a radiating `Value` condition on *every* tracer
`model_tracers` names, and `boundary_series_value` raises rather than defaulting when the prepared
file has no series for one — naming an open edge is a statement that the boundary data supplies the
whole exterior state there. So `OpenLateralBoundaryFromData` below is `oslofjorden.jl`'s, unchanged,
and this is what makes it work.

`boundary_variable_names` is extended as well as `prepare_boundary_data`: the map is what
`boundary_edge_series` reads the file back through, so a variable missing from it would be written
and then ignored. The four new entries are keyed by themselves rather than by a NorKyst variable —
nothing in the download corresponds to them, and nothing looks for them there, because
`prepare_boundary_data` only ever considers the source names listed in `parameters`.
"""
struct NPZDBoundaries{B} <: AbstractBoundaryDataConfig
    base::B
    data_root::String
    output_file::String
    plot_file::String
    open_edges::Vector{Symbol}
end

NPZDBoundaries(base) =
    NPZDBoundaries(base, base.data_root, base.output_file, base.plot_file, base.open_edges)

BoundaryData.boundary_variable_names(config::NPZDBoundaries) = merge(
    BoundaryData.boundary_variable_names(config.base),
    Dict(String(name) => String(name) for name in keys(NPZD_PROFILES)),
)

BoundaryData.download_boundary_data(target_grid, config::NPZDBoundaries) =
    download_boundary_data(target_grid, config.base)

"""
    prepare_boundary_data(target_grid, config::NPZDBoundaries; coverage)

Regrid `base`'s exterior state, then append the biogeochemical variables to the file it wrote.

Appended after the pipeline rather than inside it: `boundary_source_slab` can *derive* a variable
from others in the download, but `prepared_boundary_variable` still validates the source name
against the downloaded file, so a variable with no source at all has no place there.
"""
function BoundaryData.prepare_boundary_data(target_grid, config::NPZDBoundaries; coverage = nothing)
    result = prepare_boundary_data(target_grid, config.base; coverage)
    edges = config.open_edges
    append_boundary_variables!(result.output_file, edges, NPZD_PROFILES)

    appended = [
        BoundaryData.boundary_variable_name(edge, String(name))
        for edge in edges for name in keys(NPZD_PROFILES)
    ]

    return (; result..., variables = [result.variables; appended])
end

# --- The config itself --------------------------------------------------------------------------

# `oslofjorden.jl` is the whole science of this run — grid, bathymetry, closure, rivers, boundary
# data, atmosphere, time stepping. Everything below either passes one of its pieces through
# untouched or wraps exactly one of them. Built from the file beside this one.
base = include("oslofjorden.jl")

# `output_file` is a new name so `oslofjorden.jl`'s own `rivers_nve.nc` survives for comparison, and
# so does every other file this example writes.
rivers = NPZDRivers(
    base.river_config,
    output_file = "rivers_npzd.nc",
    plot_file = "rivers_npzd.png",
)

boundaries = NPZDBoundaries(
    NorKystBoundaryDataConfig(
        data_root = base.boundary_data_config.data_root,
        output_directory = base.boundary_data_config.output_directory,
        output_file = "boundaries_npzd.nc",
        plot_file = "boundaries_npzd.png",
        open_edges = base.boundary_data_config.open_edges,
        margin = base.boundary_data_config.margin,
        architecture = base.boundary_data_config.architecture,
        parameters = base.boundary_data_config.parameters,
        years = base.boundary_data_config.years,
    ),
)

simulation = base.simulation_config
base_model = base.model_config

FjordConfig(
    grid_config = base.grid_config,
    bathymetry_config = base.bathymetry_config,
    river_config = rivers,
    boundary_data_config = boundaries,
    atmosphere_config = base.atmosphere_config,
    # `oslofjorden.jl`'s model with the four tracers added.
    model_config = CoupledHydrostaticModel(
        buoyancy = base_model.buoyancy,
        closure = base_model.closure,
        # `oslofjorden.jl`'s scalar `WENO()`, which was already chosen so that a tracer the setup
        # does not enumerate — CATKE's `e` then, these four now — does not silently fall back to
        # an unbounded centered scheme.
        tracer_advection = base_model.tracer_advection,
        momentum_advection = base_model.momentum_advection,
        # Spelled out in full. Oceananigans merges a biogeochemical model's required tracers into
        # the model's own only when *none* of them is named already (`has_biogeochemical_tracers`),
        # and NPZD requires `:T` — its growth rate carries a Q₁₀ — which every FjordSim setup
        # names. All-but-one present is an `ArgumentError`, so the list must be complete. It is
        # also what `model_tracers` reports, which decides the forcing terms, the open boundary
        # conditions and the initial conditions, all resolved before the model exists.
        tracers = (:T, :S, :N, :P, :Z, :D),
        coriolis = base_model.coriolis,
        sea_ice = base_model.sea_ice,
        free_surface = base_model.free_surface,
        extra_kwargs = base_model.extra_kwargs,
    ),
    # Overloads model_biogeochemistry, which build_simulation calls once the grid exists.
    biogeochemistry_config = NPZDBiogeochemistry(surface_PAR = surface_par, scale_negatives = true),
    # Unchanged. The air-sea flux piece gives a top boundary condition to u, v, T and S only, so
    # the four biogeochemical tracers get Oceananigans' default no-flux surface — right for a
    # first NPZD, which has no gas exchange and no atmospheric nitrogen deposition. The open
    # southern edge is what `NPZDBoundaries` exists to feed.
    boundary_condition_configs = base.boundary_condition_configs,
    simulation_config = SimulationConfig(
        # Its own root. Checkpoints are shared per `results_root`, so pointing this at
        # `oslofjorden.jl`'s would have the two runs pick each other's up.
        results_root = fjord_results_root("oslofjorden_npzd"),
        architecture = simulation.architecture,
        writers = (
            SnapshotWriter(
                name = :ocean,
                output_file = "snapshots_ocean.nc",
                # The four new tracers beside `oslofjorden.jl`'s five. `PAR` is a biogeochemical
                # *auxiliary* field rather than a model field, so it is not nameable here.
                variables = (:T, :S, :u, :v, :e, :N, :P, :Z, :D),
                interval = 3hour,
                overwrite_existing = true,
            ),
            FieldSnapshotWriter(
                name = :surface,
                output_file = "snapshots_surface.jld2",
                variables = (:η,),
                interval = 3hour,
                overwrite_existing = true,
            ),
            CheckpointWriter(interval = 12hours, cleanup = true),
        ),
        callbacks = simulation.callbacks,
        time_stepping = simulation.time_stepping,
        # A uniform T and S, plus N, P, Z and D as their profiles on `start_date`. Not
        # `oslofjorden.jl`'s `FromInitialState()`, which reads a file and does not combine with
        # the profiles here.
        initial_conditions = merge(
            (T = 5.0, S = 33.0),
            map(profile -> (λ, φ, z) -> profile(z, simulation.start_date), NPZD_PROFILES),
        ),
        start_date = simulation.start_date,
        stop_time = simulation.stop_time,
        loops = simulation.loops,
        pickup = simulation.pickup,
    ),
)
