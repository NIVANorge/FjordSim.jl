module Plotting

export plot_bathymetry, plot_forcing, plot_rivers, plot_boundaries, plot_atmosphere, plot_path

using CairoMakie
using NCDatasets
using Oceananigans
using Oceananigans.Architectures: on_architecture
using Oceananigans.Fields: interior
using Oceananigans.Grids: x_domain, y_domain, inactive_node

using ..Configs:
    AbstractBathymetryConfig,
    AbstractForcingConfig,
    AbstractRiverConfig,
    AbstractBoundaryDataConfig,
    AbstractAtmosphereConfig
using ..Grids: open_edges

# Logical figure pixels per grid cell, and fixed margin for title/labels/colorbar. Sized from
# the grid so every cell is rendered as a distinct, unambiguous block rather than being blurred
# together by too few output pixels per cell.
const BATHYMETRY_PLOT_PIXELS_PER_CELL = 6
const BATHYMETRY_PLOT_MARGIN = (300, 150)
const FORCING_PLOT_PIXELS_PER_CELL = 3
const FORCING_PLOT_MARGIN = (200, 150)
# A boundary section is one row of cells by a handful of levels, so it gets a fixed panel size
# rather than one derived from the grid: a per-cell height would make eighteen levels unreadable.
const BOUNDARY_PLOT_PANEL_SIZE = (420, 220)
const ATMOSPHERE_PLOT_PIXELS_PER_CELL = 6
const ATMOSPHERE_PLOT_MARGIN = (200, 200)
const ATMOSPHERE_PLOT_COLUMNS = 4

# Surface fields worth eyeballing after preparation, with the colormap suiting each.
const FORCING_PLOT_PANELS = (
    (name = "T", label = "Temperature (ᵒC)", colormap = :thermal),
    (name = "S", label = "Salinity (g kg⁻¹)", colormap = :haline),
    (name = "u", label = "u (m s⁻¹)", colormap = :balance),
    (name = "v", label = "v (m s⁻¹)", colormap = :balance),
    (name = "T_lambda", label = "T relaxation rate (s⁻¹)", colormap = :viridis),
)

# Prepared boundary variables worth eyeballing, in the order they are laid out. `T`, `S`, `u` and
# `v` are sections (along the edge by depth); `eta`, `ubar` and `vbar` are lines along the edge.
const BOUNDARY_PLOT_PANELS = (
    (name = "T", label = "Temperature (ᵒC)", colormap = :thermal),
    (name = "S", label = "Salinity (g kg⁻¹)", colormap = :haline),
    (name = "u", label = "u (m s⁻¹)", colormap = :balance),
    (name = "v", label = "v (m s⁻¹)", colormap = :balance),
    (name = "eta", label = "Elevation (m)", colormap = :balance),
    (name = "ubar", label = "ubar (m s⁻¹)", colormap = :balance),
    (name = "vbar", label = "vbar (m s⁻¹)", colormap = :balance),
)

# Size of each of the two line panels beside a rivers map — discharge and temperature — since
# neither is sized per grid cell.
const RIVER_PLOT_DETAIL_SIZE = (560, 320)

# The river variable a rivers plot draws over time beside the discharge. Salinity is not worth a
# panel: a river is fresh, and every built-in source writes it as a constant zero.
const RIVER_PLOT_VALUE_NAME = "T"

# Every river gets a distinguishable colour in the two line panels. `tab20` rather than the seven
# `wong_colors` this used to be: a setup that discovers its outlets from NVE has as many rivers as
# the domain has, which is 21 on Oslofjorden, so a seven-colour cycle repeated three times. The
# cycle still repeats beyond twenty, which is why each line is also labelled.
const RIVER_PLOT_COLORS = Makie.to_colormap(:tab20)

# Every variable of a prepared atmosphere file, with the colormap suiting each.
const ATMOSPHERE_PLOT_PANELS = (
    (name = "air_temperature_2m", label = "Air temperature (K)", colormap = :thermal),
    (name = "specific_humidity_2m", label = "Specific humidity (kg kg⁻¹)", colormap = :viridis),
    (name = "air_pressure_at_sea_level", label = "Sea level pressure (Pa)", colormap = :viridis),
    (name = "precipitation", label = "Precipitation (kg m⁻² s⁻¹)", colormap = :dense),
    (name = "u_wind_10m", label = "u wind (m s⁻¹)", colormap = :balance),
    (name = "v_wind_10m", label = "v wind (m s⁻¹)", colormap = :balance),
    (name = "swrad", label = "Downwelling shortwave (W m⁻²)", colormap = :solar),
    (name = "lwrad", label = "Downwelling longwave (W m⁻²)", colormap = :solar),
)

"""
    default_figure_size(grid, pixels_per_cell, margin; columns = 1)

Figure size in logical pixels for a `columns`-panel plot of `grid`, one block per grid cell.
"""
function default_figure_size(grid, pixels_per_cell, margin; columns = 1)
    Nx, Ny, _ = size(grid)
    return (
        columns * (Nx * pixels_per_cell + margin[1]),
        Ny * pixels_per_cell + margin[2],
    )
end

"""
    plot_axes(grid, Nx, Ny)

Longitude and latitude cell coordinates spanning the domain of `grid`, for data of shape
`(Nx, Ny)`. The shape is passed in rather than taken from `grid` because a staggered forcing
variable carries one extra face row.
"""
function plot_axes(grid, Nx, Ny)
    longitude = collect(range(x_domain(grid)[1], x_domain(grid)[2], length = Nx))
    latitude = collect(range(y_domain(grid)[1], y_domain(grid)[2], length = Ny))
    return longitude, latitude
end

"""
    plot_path(config)

Where a pipeline writes its diagnostic plot. Declared here with no method: each pipeline module adds
one for its config supertype, resolving `config.plot_file` against `config.data_root`, so a new
source inherits it. A `plot_file` holding an absolute path is returned unchanged.
"""
function plot_path end

prepare_plot_file(config) = let file = plot_path(config)
    isdir(dirname(file)) || mkpath(dirname(file))
    file
end

"""
    dataset_label(config)

The dataset a config names, for plot titles: its type name without the `Config` suffix, so
`NorKystConfig` reads as "NorKyst".
"""
dataset_label(config) = replace(string(nameof(typeof(config))), r"Config$" => "")

"""
    plot_bathymetry(grid, bottom_height, config::AbstractBathymetryConfig; title, figure_size)

Write a diagnostic bathymetry plot to `plot_path(config)`: bottom height with land blanked out,
depth contours every 25 m and the coastline at `h = 0`. Returns the plot path.
"""
function plot_bathymetry(
    grid,
    bottom_height,
    config::AbstractBathymetryConfig;
    title = "Bathymetry",
    figure_size = default_figure_size(grid, BATHYMETRY_PLOT_PIXELS_PER_CELL, BATHYMETRY_PLOT_MARGIN),
)
    plot_file = prepare_plot_file(config)

    cpu_bottom_height = on_architecture(CPU(), bottom_height)
    bathymetry_data = Array(interior(cpu_bottom_height, :, :, 1))

    Nx, Ny, _ = size(grid)
    longitude, latitude = plot_axes(grid, Nx, Ny)

    figure = Figure(size = figure_size)
    axis = Axis(figure[1, 1]; xlabel = "Longitude", ylabel = "Latitude", title)

    plot = heatmap!(axis, longitude, latitude, bathymetry_data; colormap = :deep, colorrange = extrema(bathymetry_data))
    land_mask = ifelse.(bathymetry_data .>= 0, 1.0f0, NaN32)
    heatmap!(
        axis,
        longitude,
        latitude,
        land_mask;
        colormap = [:ivory, :ivory],
        colorrange = (0, 1),
        nan_color = RGBAf(0, 0, 0, 0),
    )
    Colorbar(figure[1, 2], plot; label = "Bottom height (m)")
    contour!(
        axis,
        longitude,
        latitude,
        bathymetry_data;
        levels = -collect(25:25:300),
        color = (:white, 0.35),
        linewidth = 1,
    )
    contour!(
        axis,
        longitude,
        latitude,
        bathymetry_data;
        levels = [0.0],
        color = :black,
        linewidth = 4,
    )
    save(plot_file, figure)

    return plot_file
end

"""
    surface_slice(ds, name)

The top vertical level of `name` at the first time step, as a `Float32` array with missing
values turned into `NaN` so Makie leaves land blank.
"""
function surface_slice(ds, name)
    Nz = ds.dim["Nz"]
    return Float32.(coalesce.(ds[name][:, :, Nz, 1], NaN32))
end

"""
    plot_forcing(grid, config::AbstractForcingConfig, forcing_file)

Write a diagnostic plot of the prepared forcing file at `forcing_file` to
`plot_path(config)`: one panel per `FORCING_PLOT_PANELS` entry present in the file, at the
surface level of the first time step. Returns the plot path.
"""
function plot_forcing(grid, config::AbstractForcingConfig, forcing_file)
    plot_file = prepare_plot_file(config)

    NCDataset(forcing_file) do ds
        panels = [panel for panel in FORCING_PLOT_PANELS if haskey(ds, panel.name)]
        figure = Figure(
            size = default_figure_size(
                grid,
                FORCING_PLOT_PIXELS_PER_CELL,
                FORCING_PLOT_MARGIN;
                columns = length(panels),
            ),
        )
        date = ds["time"][1]

        for (column, panel) in enumerate(panels)
            data = surface_slice(ds, panel.name)
            Nx, Ny = size(data)
            longitude, latitude = plot_axes(grid, Nx, Ny)

            axis = Axis(
                figure[1, column];
                xlabel = "Longitude",
                ylabel = column == 1 ? "Latitude" : "",
                title = panel.name,
            )
            finite = filter(isfinite, data)
            colorrange = isempty(finite) || allequal(finite) ? (0, 1) : extrema(finite)
            plot = heatmap!(
                axis,
                longitude,
                latitude,
                data;
                colormap = panel.colormap,
                colorrange,
                nan_color = :ivory,
            )
            Colorbar(figure[2, column], plot; label = panel.label, vertical = false)
        end

        Label(figure[0, :], "Prepared $(dataset_label(config)) forcing, surface level at $date", fontsize = 20)
        save(plot_file, figure)
    end

    return plot_file
end

"""
    plot_rivers(grid, config::AbstractRiverConfig, river_file, mouths)

Write a diagnostic plot of the river file at `river_file` to `plot_path(config)`.
Returns the plot path.

Three panels, one per question a river file has to answer:

- **where the water enters** — a map of the surface water with every river's columns drawn in its
  own colour and its mouth marked. The columns are the content: they are what the spreading by
  `maximum_inflow_velocity` decided, and a river that ended up across a sound or up a side arm is
  visible here and nowhere else.
- **how much** — each river's discharge over the whole axis, which is what says whether the series
  arrived and whether its floods are plausible.
- **at what temperature** — each river's `T` at its mouth over the whole axis. A river whose `T` is
  `NaN` throughout carries the temperature of the water it enters; it is named in the legend and
  draws nothing.

`mouths` comes from `add_rivers`' own result, in the order the file numbers its rivers.
"""
function plot_rivers(grid, config::AbstractRiverConfig, river_file, mouths)
    plot_file = prepare_plot_file(config)

    Nx, Ny, Nz = size(grid)
    longitude, latitude = plot_axes(grid, Nx, Ny)
    cell_width = step(range(first(longitude), last(longitude), length = Nx))
    cell_height = step(range(first(latitude), last(latitude), length = Ny))
    water = [!inactive_node(i, j, Nz, grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
    river_color(n) = RIVER_PLOT_COLORS[mod1(n, length(RIVER_PLOT_COLORS))]
    label(mouth) = "$(mouth.location.id) $(mouth.location.name)"

    NCDataset(river_file) do ds
        dates = ds["time"][:]
        discharge = Array(ds["discharge"][:, :])
        width, height = default_figure_size(grid, FORCING_PLOT_PIXELS_PER_CELL, FORCING_PLOT_MARGIN)
        figure = Figure(
            size = (
                width + RIVER_PLOT_DETAIL_SIZE[1],
                max(height, 2 * RIVER_PLOT_DETAIL_SIZE[2]),
            ),
        )

        map_axis = Axis(
            figure[1:2, 1];
            xlabel = "Longitude",
            ylabel = "Latitude",
            title = "Columns each river enters through",
        )
        heatmap!(
            map_axis, longitude, latitude, ifelse.(water, 1.0f0, NaN32);
            colormap = [:lightsteelblue, :lightsteelblue], colorrange = (0, 1), nan_color = :ivory,
        )
        for (n, mouth) in enumerate(mouths)
            cells = [
                Rect2f(
                    longitude[i] - cell_width / 2, latitude[j] - cell_height / 2,
                    cell_width, cell_height,
                ) for (i, j) in mouth.columns
            ]
            poly!(map_axis, cells; color = river_color(n))
        end
        # The mouth itself, on the cell the outlet snapped to rather than where the dataset put it —
        # the two differ by `mouth.distance`, and it is the snapped cell the columns grow from.
        scatter!(
            map_axis,
            [longitude[mouth.i] for mouth in mouths],
            [latitude[mouth.j] for mouth in mouths];
            color = [river_color(n) for n in eachindex(mouths)],
            markersize = 11,
            strokecolor = :black,
            strokewidth = 1,
        )

        discharge_axis = Axis(
            figure[1, 2];
            xlabel = "Date",
            ylabel = "Discharge (m³ s⁻¹)",
            title = "River discharge",
            yscale = Makie.pseudolog10,
        )
        for (n, mouth) in enumerate(mouths)
            lines!(discharge_axis, dates, discharge[n, :]; color = river_color(n), label = label(mouth))
        end
        # Two banks, because a discovered river list is long enough that one column of labels runs
        # off the bottom of the panel it sits in.
        axislegend(discharge_axis; position = :lt, labelsize = 8, nbanks = 2, framevisible = false)

        if haskey(ds, RIVER_PLOT_VALUE_NAME)
            value_axis = Axis(
                figure[2, 2];
                xlabel = "Date",
                ylabel = "$RIVER_PLOT_VALUE_NAME at the mouth",
                title = "River $RIVER_PLOT_VALUE_NAME",
            )
            for (n, mouth) in enumerate(mouths)
                i, j = first(mouth.columns)
                series = Float32.(coalesce.(ds[RIVER_PLOT_VALUE_NAME][i, j, :], NaN32))
                lines!(value_axis, dates, series; color = river_color(n), label = label(mouth))
            end
            axislegend(value_axis; position = :lt, labelsize = 8, nbanks = 2, framevisible = false)
        end

        Label(
            figure[0, :],
            "River inflow from $(dataset_label(config)): $(length(mouths)) river(s) over " *
            "$(sum(mouth -> length(mouth.columns), mouths)) column(s)",
            fontsize = 20,
        )
        save(plot_file, figure)
    end

    return plot_file
end

"""
    plot_boundaries(config::AbstractBoundaryDataConfig, boundary_file)

Write a diagnostic plot of the prepared boundary file at `boundary_file` to
`plot_path(config)`: one panel per `BOUNDARY_PLOT_PANELS` entry present in the file, for each edge in
the config's own `open_edges`, stacked one block of rows per edge. Returns the plot path.

Two panel shapes, because the file holds two: a full-depth variable is a section along the edge
against depth, a surface variable a line along the edge. Both are drawn against the along-edge cell
index rather than a coordinate — the file states its own dimensions but not which way round the edge
runs, and the index is what a `NaN` in a section can be traced back to.

The variables are named for their side, so the plot has to know which boundaries it is looking at; it
reads them from the config, which is where the edges are stated.
"""
function plot_boundaries(config::AbstractBoundaryDataConfig, boundary_file)
    edges = open_edges(config)
    plot_file = prepare_plot_file(config)

    NCDataset(boundary_file) do ds
        # One block of rows per open edge, so a domain open all round gets four stacked sections of the
        # same panels rather than four files.
        blocks = [(; edge, panels = boundary_panels(ds, edge)) for edge in edges]
        empty_edges = [block.edge for block in blocks if isempty(block.panels)]
        isempty(empty_edges) || error(
            "Boundary file $boundary_file carries no $(join(empty_edges, ", ")) variables to plot.",
        )

        columns = maximum(block -> length(block.panels), blocks)
        figure = Figure(
            size = (
                columns * BOUNDARY_PLOT_PANEL_SIZE[1],
                length(blocks) * BOUNDARY_PLOT_PANEL_SIZE[2] + FORCING_PLOT_MARGIN[2],
            ),
        )
        date = ds["time"][1]
        depth = haskey(ds, "Nz") ? Array{Float64}(ds["Nz"][:]) : Float64[]

        for (block_index, block) in enumerate(blocks)
            # Two rows per edge: the sections, then their colorbars.
            axis_row = 2 * block_index - 1

            for (column, panel) in enumerate(block.panels)
                data = boundary_slice(ds, panel.variable)
                axis = Axis(
                    figure[axis_row, column];
                    xlabel = "Cell along the :$(block.edge) boundary",
                    ylabel = column == 1 ? (ndims(data) == 2 ? "Depth (m)" : "") : "",
                    title = panel.name,
                )

                if ndims(data) == 2
                    finite = filter(isfinite, data)
                    colorrange = isempty(finite) || allequal(finite) ? (0, 1) : extrema(finite)
                    plot = heatmap!(
                        axis,
                        1:size(data, 1),
                        depth,
                        data;
                        colormap = panel.colormap,
                        colorrange,
                        nan_color = :ivory,
                    )
                    Colorbar(figure[axis_row+1, column], plot; label = panel.label, vertical = false)
                else
                    lines!(axis, 1:length(data), data)
                    axis.ylabel = panel.label
                end
            end
        end

        Label(
            figure[0, :],
            "Prepared $(dataset_label(config)) $(join((string(":", edge) for edge in edges), ", ")) " *
            "boundary at $date",
            fontsize = 20,
        )
        save(plot_file, figure)
    end

    return plot_file
end

"""
    boundary_panels(ds, edge)

The `BOUNDARY_PLOT_PANELS` entries the file carries for one edge, each with the prefixed variable name
it is read from.
"""
boundary_panels(ds, edge) = [
    (; panel..., variable = string(edge, "_", panel.name))
    for panel in BOUNDARY_PLOT_PANELS if haskey(ds, string(edge, "_", panel.name))
]

"""
    boundary_slice(ds, name)

A prepared boundary variable at the first time step, as a `Float32` array with missing values turned
into `NaN` so Makie leaves land blank: a `(along, Nz)` section for a full-depth variable, a
`(along,)` line for a surface one.
"""
function boundary_slice(ds, name)
    variable = ds[name]
    slice = variable[ntuple(_ -> Colon(), ndims(variable) - 1)..., 1]
    return Float32.(coalesce.(slice, NaN32))
end

"""
    plot_atmosphere(config::AbstractAtmosphereConfig, atmosphere_file)

Write a diagnostic plot of the prepared atmosphere file at `atmosphere_file` to
`plot_path(config)`: one panel per `ATMOSPHERE_PLOT_PANELS` entry present in the file, at the first
time step. Returns the plot path.

Unlike `plot_forcing` this takes no grid, because the prepared file carries its own `lon` and `lat`
axes — the atmosphere grid is independent of the ocean grid.
"""
function plot_atmosphere(config::AbstractAtmosphereConfig, atmosphere_file)
    plot_file = prepare_plot_file(config)

    NCDataset(atmosphere_file) do ds
        panels = [panel for panel in ATMOSPHERE_PLOT_PANELS if haskey(ds, panel.name)]
        isempty(panels) && error("$atmosphere_file holds none of the expected atmosphere variables")

        longitude = Array(ds["lon"][:])
        latitude = Array(ds["lat"][:])
        columns = min(length(panels), ATMOSPHERE_PLOT_COLUMNS)
        rows = cld(length(panels), columns)
        figure = Figure(
            size = (
                columns * (length(longitude) * ATMOSPHERE_PLOT_PIXELS_PER_CELL + ATMOSPHERE_PLOT_MARGIN[1]),
                rows * (length(latitude) * ATMOSPHERE_PLOT_PIXELS_PER_CELL + ATMOSPHERE_PLOT_MARGIN[2]),
            ),
        )
        date = ds["time"][1]

        for (position, panel) in enumerate(panels)
            # Two figure rows per panel row: the heatmap and its horizontal colorbar.
            row = 2 * cld(position, columns) - 1
            column = mod1(position, columns)
            data = Float32.(coalesce.(ds[panel.name][:, :, 1], NaN32))

            axis = Axis(
                figure[row, column];
                xlabel = "Longitude",
                ylabel = column == 1 ? "Latitude" : "",
                title = panel.name,
            )
            finite = filter(isfinite, data)
            colorrange = isempty(finite) || allequal(finite) ? (0, 1) : extrema(finite)
            plot = heatmap!(
                axis,
                longitude,
                latitude,
                data;
                colormap = panel.colormap,
                colorrange,
                nan_color = :ivory,
            )
            Colorbar(figure[row+1, column], plot; label = panel.label, vertical = false)
        end

        Label(figure[0, :], "Prepared $(dataset_label(config)) atmosphere at $date", fontsize = 20)
        save(plot_file, figure)
    end

    return plot_file
end

end  # module Plotting
