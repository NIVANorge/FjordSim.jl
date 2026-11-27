# The Oslofjord setup, 160 m, October 2024 to October 2026. What it covers, how to prepare and run
# it, and why its values are what they are: README.md beside this file.

using FjordSim
using Oceananigans:
    CATKEVerticalDiffusivity,
    HydrostaticSphericalCoriolis,
    SeawaterBuoyancy,
    WENO,
    WENOVectorInvariant
using Oceananigans.TurbulenceClosures: HorizontalScalarBiharmonicDiffusivity
using Oceananigans.Units: day, days, hour, hours, minutes, second
using Dates: DateTime, Day, Month, Second
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using FjordSim.Observations: VannmiljoObservations
using FjordSim.Validation: ValidationConfig
using NumericalEarth: FreezingLimitedOceanTemperature

include("oslofjorden_hooks.jl")

# --- Stations ---------------------------------------------------------------------------------
#
# Each group states where its positions come from.

# Kartverket's permanent tide gauges, from its station register. Exact.
const OSLOFJORD_TIDE_GAUGES = [
    Station(name = "Viker",      longitude = 10.949769, latitude = 59.036046),
    Station(name = "Oscarsborg", longitude = 10.604861, latitude = 59.678073),
    Station(name = "Oslo",       longitude = 10.734510, latitude = 59.908559),
]

# NIVA's Ytre Oslofjord stations, METreport 11/2017 Table 2. Exact. AquaMonitor puts OF-1 5 km west,
# in 76 m of water, where its casts reach 440 m; the table's position is 462 m deep.
const OSLOFJORD_CTD_STATIONS = [
    Station(name = "D-2",  longitude = 10.4210, latitude = 59.6280),  # Inner Drammensfjord
    Station(name = "D-3",  longitude = 10.3140, latitude = 59.7060),  # Solumstrand
    Station(name = "LA-1", longitude = 10.0520, latitude = 59.0190),  # Larviksfjord
    Station(name = "MO-2", longitude = 10.6780, latitude = 59.4840),  # Kippenes
    Station(name = "OF-1", longitude = 10.7540, latitude = 59.0410),  # Torbjørnskjær
    Station(name = "OF-5", longitude = 10.4580, latitude = 59.4870),  # Breiangen
    Station(name = "S-9",  longitude = 11.1620, latitude = 59.1140),  # Haslau, Singlefjord
    Station(name = "TØ-1", longitude = 10.3550, latitude = 59.2030),  # Vestfjord
    Station(name = "Ø-1",  longitude = 10.8340, latitude = 59.1370),  # Leira, Vesterelva
]

# The programme's other stations with casts in 2024-2026, at AquaMonitor's positions. Exact.
const OSLOFJORD_PROGRAMME_STATIONS = [
    Station(name = "BO-1",  longitude = 10.536233, latitude = 59.189933),  # Bolærne,          cast to 58 m
    Station(name = "I-1",   longitude = 11.001966, latitude = 59.109432),  # Ramsø, Østerelva, cast to 51 m
    Station(name = "I-4",   longitude = 10.951380, latitude = 59.184500),  # Kallera,          cast to 11 m
    Station(name = "KF-1",  longitude = 10.729045, latitude = 59.288682),  # Krokstadfjorden,  cast to 15 m
    Station(name = "SAN-5", longitude = 10.269300, latitude = 59.549400),  # Sandebukta,       cast to 51 m
    Station(name = "SF-3",  longitude = 10.247500, latitude = 59.068000),  # Kvernberget,      cast to 83 m
    Station(name = "VAL-1", longitude = 10.525400, latitude = 59.270100),  # Vallø,            cast to 129 m
]

# NIVA's Inner Oslofjord stations, at AquaMonitor's positions. Exact.
const INNER_OSLOFJORD_STATIONS = [
    Station(name = "Dk1", longitude = 10.569384, latitude = 59.814999),  # Steilene,         cast to 99 m
    Station(name = "Ep1", longitude = 10.723783, latitude = 59.786301),  # Bunnefjorden,     cast to 149 m
    Station(name = "Im2", longitude = 10.628217, latitude = 59.622032),  # Elle,             cast to 202 m
    Station(name = "Fl1", longitude = 10.574734, latitude = 59.754066),  # Vestfjorden,      cast to 166 m
    Station(name = "Gk1", longitude = 10.543083, latitude = 59.701218),  # Gråøyrenna,       cast to 112 m
    Station(name = "Bl4", longitude = 10.568850, latitude = 59.880833),  # Bærumsbassenget,  cast to 31 m
    Station(name = "Cj1", longitude = 10.506983, latitude = 59.840483),  # Holmenfjorden,    cast to 56 m
]

# The Drammensfjord's deepest path, from the Drammenselva mouth over the Svelvik sill to Breiangen
# (README.md): distance along it and depth on the 193 m grid it was found on.
const DRAMMENSFJORD_SECTION = [
    Station(name = "DS-00", longitude = 10.21026, latitude = 59.73954),  #  0.00 km,   5.0 m
    Station(name = "DS-01", longitude = 10.21709, latitude = 59.73597),  #  0.55 km,   4.3 m
    Station(name = "DS-02", longitude = 10.22393, latitude = 59.73418),  #  1.02 km,  13.0 m
    Station(name = "DS-03", longitude = 10.23077, latitude = 59.73239),  #  1.49 km,   6.2 m
    Station(name = "DS-04", longitude = 10.23761, latitude = 59.73239),  #  1.88 km,   7.0 m
    Station(name = "DS-05", longitude = 10.24444, latitude = 59.73418),  #  2.35 km,  35.7 m
    Station(name = "DS-06", longitude = 10.25128, latitude = 59.73418),  #  2.73 km,  44.1 m
    Station(name = "DS-07", longitude = 10.25812, latitude = 59.73239),  #  3.20 km,  52.8 m
    Station(name = "DS-08", longitude = 10.26496, latitude = 59.72881),  #  3.76 km,  64.0 m
    Station(name = "DS-09", longitude = 10.27179, latitude = 59.72524),  #  4.31 km,  70.1 m
    Station(name = "DS-10", longitude = 10.27863, latitude = 59.72166),  #  4.86 km,  76.4 m
    Station(name = "DS-11", longitude = 10.28547, latitude = 59.71808),  #  5.42 km,  79.6 m
    Station(name = "DS-12", longitude = 10.29231, latitude = 59.71451),  #  5.97 km,  83.2 m
    Station(name = "DS-13", longitude = 10.29915, latitude = 59.71093),  #  6.53 km,  86.5 m
    Station(name = "DS-14", longitude = 10.30598, latitude = 59.70735),  #  7.08 km,  94.6 m
    Station(name = "DS-15", longitude = 10.31282, latitude = 59.70735),  #  7.47 km,  96.6 m
    Station(name = "DS-16", longitude = 10.31966, latitude = 59.70557),  #  7.94 km,  98.4 m
    Station(name = "DS-17", longitude = 10.32650, latitude = 59.70557),  #  8.32 km,  99.4 m
    Station(name = "DS-18", longitude = 10.33333, latitude = 59.70378),  #  8.79 km, 102.8 m
    Station(name = "DS-19", longitude = 10.34017, latitude = 59.70378),  #  9.18 km,  95.6 m
    Station(name = "DS-20", longitude = 10.34701, latitude = 59.70378),  #  9.56 km,  64.4 m
    Station(name = "DS-21", longitude = 10.35385, latitude = 59.70378),  #  9.95 km,  90.8 m
    Station(name = "DS-22", longitude = 10.36068, latitude = 59.70020),  # 10.50 km, 108.8 m
    Station(name = "DS-23", longitude = 10.36752, latitude = 59.69662),  # 11.05 km, 110.9 m
    Station(name = "DS-24", longitude = 10.37436, latitude = 59.69305),  # 11.61 km, 113.4 m
    Station(name = "DS-25", longitude = 10.38120, latitude = 59.68947),  # 12.16 km, 114.3 m
    Station(name = "DS-26", longitude = 10.38803, latitude = 59.68589),  # 12.72 km, 114.1 m
    Station(name = "DS-27", longitude = 10.38803, latitude = 59.68232),  # 13.11 km, 114.0 m
    Station(name = "DS-28", longitude = 10.39145, latitude = 59.67874),  # 13.59 km, 114.0 m
    Station(name = "DS-29", longitude = 10.39829, latitude = 59.67516),  # 14.14 km, 114.7 m
    Station(name = "DS-30", longitude = 10.39829, latitude = 59.67159),  # 14.54 km, 114.0 m
    Station(name = "DS-31", longitude = 10.40171, latitude = 59.66801),  # 15.02 km, 114.0 m
    Station(name = "DS-32", longitude = 10.40855, latitude = 59.66443),  # 15.57 km, 118.3 m
    Station(name = "DS-33", longitude = 10.41197, latitude = 59.66086),  # 16.05 km, 119.9 m
    Station(name = "DS-34", longitude = 10.41197, latitude = 59.65728),  # 16.45 km, 120.8 m
    Station(name = "DS-35", longitude = 10.41197, latitude = 59.65370),  # 16.84 km, 121.1 m
    Station(name = "DS-36", longitude = 10.41197, latitude = 59.65013),  # 17.24 km, 121.6 m
    Station(name = "DS-37", longitude = 10.41197, latitude = 59.64655),  # 17.64 km, 121.9 m
    Station(name = "DS-38", longitude = 10.41197, latitude = 59.64297),  # 18.04 km, 121.9 m
    Station(name = "DS-39", longitude = 10.41538, latitude = 59.63940),  # 18.51 km, 121.9 m
    Station(name = "DS-40", longitude = 10.41880, latitude = 59.63582),  # 18.99 km, 121.9 m
    Station(name = "DS-41", longitude = 10.42222, latitude = 59.63224),  # 19.47 km, 116.8 m
    Station(name = "DS-42", longitude = 10.42222, latitude = 59.62867),  # 19.86 km, 116.8 m
    Station(name = "DS-43", longitude = 10.42222, latitude = 59.62509),  # 20.26 km, 114.8 m
    Station(name = "DS-44", longitude = 10.41538, latitude = 59.62151),  # 20.82 km,  85.0 m
    Station(name = "DS-45", longitude = 10.41197, latitude = 59.61973),  # 21.21 km,  33.9 m
    Station(name = "DS-46", longitude = 10.41197, latitude = 59.61615),  # 21.61 km,  13.5 m
    Station(name = "DS-47", longitude = 10.41197, latitude = 59.61257),  # 22.00 km,  11.6 m
    Station(name = "DS-48", longitude = 10.41197, latitude = 59.60900),  # 22.40 km,  11.6 m
    Station(name = "DS-49", longitude = 10.41538, latitude = 59.60542),  # 22.88 km,  11.6 m
    Station(name = "DS-50", longitude = 10.41880, latitude = 59.60184),  # 23.36 km,  12.4 m
    Station(name = "DS-51", longitude = 10.42222, latitude = 59.59827),  # 23.83 km,  12.2 m
    Station(name = "DS-52", longitude = 10.42564, latitude = 59.59469),  # 24.31 km,  12.5 m
    Station(name = "DS-53", longitude = 10.42564, latitude = 59.59290),  # 24.51 km,  10.8 m, Svelvik sill
    Station(name = "DS-54", longitude = 10.42564, latitude = 59.59111),  # 24.71 km,  12.3 m
    Station(name = "DS-55", longitude = 10.42222, latitude = 59.58754),  # 25.18 km,  14.4 m
    Station(name = "DS-56", longitude = 10.42222, latitude = 59.58396),  # 25.58 km,  17.0 m
    Station(name = "DS-57", longitude = 10.42222, latitude = 59.58038),  # 25.98 km,  18.1 m
    Station(name = "DS-58", longitude = 10.42222, latitude = 59.57681),  # 26.38 km,  19.0 m
    Station(name = "DS-59", longitude = 10.42222, latitude = 59.57323),  # 26.77 km,  22.5 m
    Station(name = "DS-60", longitude = 10.41880, latitude = 59.56965),  # 27.25 km,  27.7 m
    Station(name = "DS-61", longitude = 10.41197, latitude = 59.56608),  # 27.80 km,  35.7 m
    Station(name = "DS-62", longitude = 10.41197, latitude = 59.56250),  # 28.20 km,  25.0 m
    Station(name = "DS-63", longitude = 10.41197, latitude = 59.55892),  # 28.60 km,  22.7 m
    Station(name = "DS-64", longitude = 10.41197, latitude = 59.55535),  # 29.00 km,  22.7 m
    Station(name = "DS-65", longitude = 10.41197, latitude = 59.55177),  # 29.40 km,  23.9 m
    Station(name = "DS-66", longitude = 10.41197, latitude = 59.54819),  # 29.79 km,  27.6 m
    Station(name = "DS-67", longitude = 10.40855, latitude = 59.54462),  # 30.27 km,  32.0 m
    Station(name = "DS-68", longitude = 10.40855, latitude = 59.54104),  # 30.67 km,  37.4 m
    Station(name = "DS-69", longitude = 10.40513, latitude = 59.53746),  # 31.15 km,  56.0 m
    Station(name = "DS-70", longitude = 10.40513, latitude = 59.53389),  # 31.54 km,  73.6 m
    Station(name = "DS-71", longitude = 10.40171, latitude = 59.53031),  # 32.02 km,  88.9 m
    Station(name = "DS-72", longitude = 10.40171, latitude = 59.52673),  # 32.42 km, 111.0 m
    Station(name = "DS-73", longitude = 10.40171, latitude = 59.52316),  # 32.82 km, 125.6 m
    Station(name = "DS-74", longitude = 10.40513, latitude = 59.51958),  # 33.29 km, 157.2 m
    Station(name = "DS-75", longitude = 10.40513, latitude = 59.51600),  # 33.69 km, 156.0 m
    Station(name = "DS-76", longitude = 10.40513, latitude = 59.51422),  # 33.89 km, 153.4 m
]

# IMR's ØkoKyst and Tisler stations, at Vannmiljø's positions. Exact. Torbjørnskjær, Breiangen vest,
# Bastøy and Missingene are NIVA's OF-1, OF-5, OF-4 and OF-2.
const IMR_STATIONS = [
    Station(name = "Bastøy",         longitude = 10.590569, latitude = 59.358670),  # cast to 280 m
    Station(name = "Breiangen vest", longitude = 10.458333, latitude = 59.486668),  # cast to 190 m
    Station(name = "Færder",         longitude = 10.550012, latitude = 58.966672),  # cast to 150 m
    Station(name = "Missingene",     longitude = 10.653178, latitude = 59.191948),  # cast to 350 m
    Station(name = "Tisler-1",       longitude = 10.975667, latitude = 58.992333),  # cast to 130 m
    Station(name = "Tisler-5",       longitude = 10.958000, latitude = 58.997333),  # cast to 120 m
    Station(name = "Torbjørnskjær",  longitude = 10.760812, latitude = 59.040700),  # cast to 440 m
]

const NIVA_CTD_STATIONS =
    vcat(OSLOFJORD_CTD_STATIONS, OSLOFJORD_PROGRAMME_STATIONS, INNER_OSLOFJORD_STATIONS)

# The NIVA stations whose casts start the water below the Drøbak and Svelvik sills
# (`CastCorrectedInitialState`), at AquaMonitor's positions. Exact.
const OSLOFJORD_BASIN_CAST_STATIONS = [
    Station(name = "Im2", longitude = 10.628217, latitude = 59.622032),  # Elle, Drøbak Sound
    Station(name = "D-2", longitude = 10.420983, latitude = 59.627766),  # Drammensfjord
    Station(name = "D-3", longitude = 10.313517, latitude = 59.705883),  # Drammensfjord
    Station(name = "Hm4", longitude = 10.612483, latitude = 59.675934),  # Oscarsborg
    Station(name = "En-1", longitude = 10.553669, latitude = 59.679659), # Engene
    Station(name = "Gk1", longitude = 10.543083, latitude = 59.701218),  # Gråøyrenna
    Station(name = "Gl2", longitude = 10.579634, latitude = 59.717384),  # Håøya
    Station(name = "Fl1", longitude = 10.574734, latitude = 59.754066),  # Vestfjorden
    Station(name = "Ej1", longitude = 10.516217, latitude = 59.793167),  # VEAS
    Station(name = "Dk1", longitude = 10.569384, latitude = 59.814999),  # Steilene
    Station(name = "Dm1", longitude = 10.617300, latitude = 59.827200),  # Steilene Nord
    Station(name = "Ep1", longitude = 10.723783, latitude = 59.786301),  # Bunnefjorden
    Station(name = "Cp2", longitude = 10.701384, latitude = 59.865002),  # Oksvalflua
    Station(name = "Cq1", longitude = 10.736667, latitude = 59.872849),  # Bekkelagsbassenget
    Station(name = "Bn1", longitude = 10.646733, latitude = 59.880634),  # Lysakerfjorden
]

# Viker's tide over the boundary's, from `scripts/measure_boundary_tide.jl` over 2024-10-01 to
# 2026-09-30.
const VIKER_TIDE_GAINS = [
    TidalGain("M2", 0.704, 3.3),     # 17.41 cm at the boundary, 12.26 observed, ± 0.08
    TidalGain("S2", 1.377, -27.6),   # 1.72 cm at the boundary, 2.37 observed, ± 0.08
    TidalGain("N2", 2.05, 1.7),      # 1.44 cm at the boundary, 2.95 observed, ± 0.08
    TidalGain("K2", 1.315, -10.0),   # 0.52 cm at the boundary, 0.68 observed, ± 0.08
    TidalGain("O1", 2.175, 21.0),    # 1.24 cm at the boundary, 2.7 observed, ± 0.08
    TidalGain("K1", 0.244, -95.1),   # 1.11 cm at the boundary, 0.27 observed, ± 0.08
]

data_root = fjord_data_root("oslofjorden")
start_date = DateTime(2024, 10, 1)
# Two years, two whole winters.
stop_date = DateTime(2026, 10, 1)
# The full-field writers start a new file each month, so a finished month can leave the VM.
monthly_splits = [(start_date + Month(month) - start_date) / Second(1) for month in 1:23]
FT = Float32
# A geometric stretch from a 1 m surface cell to the deepest sounding, with 41-132 m in 5.06 m levels
# against WENO's vertical diffusion in the sill basins (README.md).
z_faces = vcat(
    [-467.0, -433.5, -400.0, -366.5, -333.0, -299.5, -266.0, -232.5, -199.0, -165.5],
    range(-132.0, -41.0; length = 19),
    [-32.0, -25.0, -19.0, -14.5, -10.8, -7.9, -5.5, -3.7, -2.2, -1.0, 0.0],
)

FjordConfig(
    grid_config = EvenGrid(
        # 160 m both ways at 59.44°N. Every coefficient below tied to Δx is scaled to it.
        size      = (424, 681, length(z_faces) - 1),
        halo      = (7, 7, 7),
        longitude = (10.00, 11.20),
        latitude  = (58.95, 59.93),
        z_faces   = z_faces,
        float_type = FT,
    ),
    # Basin and sill depths each knob gives: README.md.
    bathymetry_config = DybdedataConfig(
        data_root             = data_root,
        output_file           = "bathymetry.nc",
        plot_file             = "bathymetry.png",
        # A 40 m native grid, which holds Svelvik's 11 m.
        raw_resolution_factor = 4,
        padding_cells         = 2,
        include_contours      = false,
        contour_stride        = 10,
        interpolation_passes  = 1,
        major_basins          = 1,
        minimum_depth         = 2.0,
        # No headland within five rows of the open edge.
        open_boundary_land_cells = 5,
        # Every island kept, Kaholmen in the Drøbak Sound among them.
        max_island_cells      = 0,
        close_narrow_passages = true,
        spike_ratio           = 0.5,
        minimum_cell_fraction = 0.2,
        # 0.25 took 34 m off Fl1 and 24 m off Hm4; 0.5 takes 8 and 4.
        max_slope_factor      = 0.5,
        # Regridding took Svelvik to 5.6 m against 11 m real.
        preserve_sills        = true,
        # The Drøbak jetty, whose ~6 m openings the grid cannot hold. Without it the sill is the
        # 23.5 m channel the jetty closes rather than the 21 m east channel.
        depth_caps            = [
            DepthCap(longitude = [10.6096, 10.6100, 10.6005], latitude = [59.6725, 59.6648, 59.6641], depth = 6.0),
        ],
        # Walled-off cells in the boundary row sat below −1.5 °C.
        raise_boundary_pockets = true,
        geonorge_cache        = true,
        regrid_cache          = false,
        # Removes the water the gridding put where the FileGDB has no data (README.md).
        max_sounding_distance = 2000.0,
    ),
    # The 160 m nest where it reaches, the 800 m run elsewhere, and the sill basins from the NIVA
    # casts within 14 days of the start, Im2's only below its sill (README.md).
    initial_state_config = CastCorrectedInitialState(
        NorKystV3InitialStateConfig(
            data_root        = data_root,
            output_directory = "norkyst_v3_initial",
            output_file      = "initial_state.nc",
            nest_member      = "m71",
            parent_member    = "m00",
            parameters       = ["temperature", "salinity", "u_eastward", "v_northward"],
            start_date       = start_date,
        ),
        VannmiljoObservations(
            data_root        = data_root,
            output_directory = joinpath("observations", "vannmiljo_niva_initial"),
            start_date       = start_date - Day(14),
            stop_date        = start_date + Day(14),
            contractors      = ["NIVA"],
            stations         = OSLOFJORD_BASIN_CAST_STATIONS,
        ),
        Day(14),
        ["Im2"],
    ),
    river_config = NVERiversConfig(
        data_root  = data_root,
        output_file = "rivers_nve.nc",
        plot_file   = "rivers_nve.png",
        start_date  = start_date,
        stop_date   = stop_date,
        # Every mouth this size or larger, discovered from NVE's river network (README.md).
        minimum_discharge = 0.5,
        # A vertical Courant number of 0.18 on the 1 m top cell at the 3-minute time step.
        maximum_inflow_velocity = 1.0e-3,
        # Gauges attached to discovered mouths, keyed by `vassdragsnr`. The rest run at their
        # catchment's normal.
        outlets = [
            # Glomma's two mouths share Solbergfoss's discharge, split by assumption (README.md),
            # with Sarpfossen's temperature.
            NVERiver(
                vassdragsnr = "002.A21", name = "Glomma (Osterelva)",
                discharge_station = "2.605.0", temperature_station = "2.1087.0",
                discharge_fraction = 2 // 3,
            ),
            NVERiver(
                vassdragsnr = "002.2A", name = "Glomma (Vesterelva)",
                discharge_station = "2.605.0", temperature_station = "2.1087.0",
                discharge_fraction = 1 // 3,
            ),
            NVERiver(
                vassdragsnr = "012.A2", name = "Drammenselva",
                discharge_station = "12.534.0", temperature_station = "12.534.0",
            ),
            # Mossevassdraget. No temperature series.
            NVERiver(vassdragsnr = "003.A4", discharge_station = "3.23.0"),
            # Akerselva's temperature gauge is a sun-exposed sensor, so none is attached.
            NVERiver(
                vassdragsnr = "006.A10", name = "Akerselva",
                discharge_station = "6.38.0",
            ),
            # Lysakerelva, as Akerselva.
            NVERiver(vassdragsnr = "007.A0", discharge_station = "7.29.0"),
            NVERiver(
                vassdragsnr = "008.A2",
                discharge_station = "8.2.0", temperature_station = "8.2.0",
            ),
            # Lierelva. No temperature series.
            NVERiver(vassdragsnr = "011.A0", discharge_station = "11.6.0"),
        ],
    ),
    # Hourly, so the Flather condition sees the tide. The reference run `m00`, the one member
    # published over the whole period; its levels stop at 300 m, the edge reaches 317 m.
    boundary_data_config = FerryBoxCorrectedBoundaries(
        TideCorrectedBoundaries(
            NorKystV3BoundaryDataConfig(
                data_root        = data_root,
                output_directory = "norkyst_v3_hourly",
                output_file      = "boundaries.nc",
                plot_file        = "boundaries.png",
                # The western and eastern walls stay closed (README.md).
                open_edges       = :south,
                margin           = 0.05,
                architecture     = :auto,
                member           = "m00",
                # The download writes `ubar_eastward` and `vbar_northward` as depth means.
                parameters       = [
                    "temperature", "salinity", "u_eastward", "v_northward", "zeta",
                    "ubar_eastward", "vbar_northward",
                ],
                start_date       = start_date,
                stop_date        = stop_date,
            ),
            VIKER_TIDE_GAINS,
        );
        # NIVA's FerryBox minutes (README.md).
        directory      = joinpath(data_root, "observations", "niva_ferrybox_2024-2026"),
        # Minutes within 0.05° (5.6 km) north of the edge's cell centres.
        band           = 0.05,
        # A week, so a day with one crossing does not set the edge alone.
        smoothing_days = 7,
        # Fitted to the Færder casts (README.md).
        taper_depth    = 25.0,
    ),
    atmosphere_config = NORA3Config(
        data_root        = data_root,
        output_directory = "nora3",
        output_file      = "atmosphere.nc",
        plot_file        = "atmosphere.png",
        resolution       = 0.02,
        padding          = 0.1,
        start_date       = start_date,
        stop_date        = stop_date,
        # Snow falls into the water, which supplies the heat that melts it.
        snowfall         = true,
    ),
    model_config = CoupledHydrostaticModel(
        buoyancy           = SeawaterBuoyancy(FT, equation_of_state = TEOS10EquationOfState(FT)),
        # A sponge damps what the open boundary radiates (README.md).
        closure            = BoundarySponge(
            # `minimum_tke` is the vertical mixing below the surface layer (README.md).
            base = (
                CATKEVerticalDiffusivity(FT; minimum_tke = 7e-7),
                # The 2Δx mode damped in 72 min; κ = 0, as diffusion overshot at fronts (README.md).
                HorizontalScalarBiharmonicDiffusivity(FT; ν = 9.4e3, κ = 0),
            ),
            width_cells = 20,
            # 30 and 15 on 193 m cells, scaled by (160/193)².
            viscosity   = 20.0,
            diffusivity = 10.0,
        ),
        # Upwind at the river mouths; T's WENO falls back to upwind near walls in the top four cells
        # only; S's WENO is bounded (README.md). `e` is named, or it gets `Centered()`.
        tracer_advection   = (
            T = RiverMouthAdvection(
                interior = SurfaceLayerAdvection(
                    surface = WENO(FT; minimum_buffer_upwind_order = 1),
                    interior = WENO(FT),
                    surface_levels = 4,
                ),
                width_cells = 3,
                narrow_cells = false,
            ),
            S = RiverMouthAdvection(
                interior = WENO(FT; bounds = (0, 40)),
                width_cells = 3,
                narrow_cells = false,
            ),
            e = WENO(FT; minimum_buffer_upwind_order = 1),
        ),
        momentum_advection = WENOVectorInvariant(FT),
        tracers            = (:T, :S),
        coriolis           = HydrostaticSphericalCoriolis(FT),
        # Clamped at freezing, with no ice (README.md).
        sea_ice            = FreezingLimitedOceanTemperature(FT),
        # The blue-green band's 3 m from the fjord's 4-6 m Secchi depths (README.md).
        radiative_forcing  = TwoColorRadiationConfig(
            first_color_fraction = 0.58,
            first_decay_scale    = 0.35,
            second_decay_scale   = 3.0,
        ),
        free_surface       = SplitExplicitFreeSurfaceConfig(cfl = 0.7),
        extra_kwargs       = (;),
    ),
    # Nudged hard on inflow, radiating on outflow (Marchesiello et al., 2001).
    boundary_condition_configs = (
        AirSeaFluxes(),
        QuadraticBottomDrag(coefficient = 0.003),
        OpenLateralBoundaryFromData(
            inflow_timescale  = 3hours,
            outflow_timescale = 360days,
        ),
    ),
    run_config = RunConfig(
        results_root       = fjord_results_root("oslofjorden"),
        architecture       = :auto,
        writers = (
            # Daily, a file a month; `e` for diagnosing the mixing.
            SnapshotWriter(
                name               = :ocean,
                output_file        = "snapshots_ocean.nc",
                variables          = (:T, :S, :u, :v, :e),
                interval           = 1day,
                split_times        = monthly_splits,
                overwrite_existing = true,
            ),
            # The 3.7-5.5 m level, where the FerryBox's intake draws. JLD2, as η cannot go in a
            # NetCDF writer (`FieldSnapshotWriter`).
            SurfaceSnapshotWriter(
                name               = :surface,
                output_file        = "snapshots_surface.jld2",
                variables          = (:η, :T, :S, :u, :v),
                depth              = 4,
                interval           = 1hour,
                split_times        = monthly_splits,
                overwrite_existing = true,
            ),
            # Hourly gives the shortest constituent analysed, MS4 (6.1 h), six samples a period.
            FieldStationWriter(
                name               = :tidegauge,
                output_file        = "stations_tidegauge.jld2",
                variables          = (:η,),
                stations           = OSLOFJORD_TIDE_GAUGES,
                interval           = 1hour,
                overwrite_existing = true,
            ),
            SectionWriter(
                name               = :section,
                output_file        = "stations_section.nc",
                variables          = (:T, :S, :e, :u, :v),
                stations           = DRAMMENSFJORD_SECTION,
                interval           = 1hour,
                overwrite_existing = true,
            ),
            CheckpointWriter(interval = 12hours, cleanup = true),
        ),
        # The downslope flow takes dense water over the Drøbak sill down the slope to its own level
        # rather than letting it convect (README.md).
        callbacks = (
            ProgressCallback(name = :progress, interval = 1hour, report = progress),
            DownslopeFlow(name = :downslope, tracers = (:T, :S), friction_timescale = 1000),
        ),
        time_stepping = AdaptiveTimeStep(
            initial_time_step    = 1second,
            cfl                  = 0.3,
            max_time_step        = 3minutes,
            max_time_step_change = 1.01,
        ),
        initial_conditions = FromInitialState(),
        start_date         = start_date,
        stop_date          = stop_date,
        loops              = 1,
        pickup             = false,
    ),
    # One source per contractor, so each cast is scored once. Vannmiljø places Ø-1 and Cj1 about
    # 600 m from the positions above, hence 1 km.
    observation_configs = (
        VannmiljoObservations(
            data_root        = data_root,
            output_directory = joinpath("observations", "vannmiljo_niva"),
            start_date       = start_date,
            stop_date        = stop_date,
            contractors      = ["NIVA"],
            stations         = NIVA_CTD_STATIONS,
            match_distance   = 1000.0,
            source           = "NIVA",
        ),
        VannmiljoObservations(
            data_root        = data_root,
            output_directory = joinpath("observations", "vannmiljo_imr"),
            start_date       = start_date,
            stop_date        = stop_date,
            contractors      = ["Havforskningsinstituttet"],
            stations         = IMR_STATIONS,
            source           = "IMR",
        ),
    ),
    # Below -1 °C is colder than the freezing clamp allows; above 36 g/kg saltier than any inflow.
    validation_config = ValidationConfig(
        spin_up          = Day(31),
        depth_bands      = [(0, 5), (5, 20), (20, 50), (50, Inf)],
        station_groups   = [
            "Outer fjord" => [
                s.name for s in vcat(OSLOFJORD_CTD_STATIONS, OSLOFJORD_PROGRAMME_STATIONS, IMR_STATIONS)
                if !(s.name in ("D-2", "D-3"))
            ],
            "Drammensfjord" => ["D-2", "D-3"],
            "Inner Oslofjord" => [s.name for s in INNER_OSLOFJORD_STATIONS],
        ],
        search_radius    = 3,
        plausible_ranges = Dict("T" => (-1.0, 25.0), "S" => (0.0, 36.0)),
        stratification   = ("S", 1.0, 25.0),
    ),
)
