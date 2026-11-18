# Validation figures. Every function takes series and metrics rather than configs or paths, for the
# same reason `metrics.jl` does: the figure that compares a tide gauge is the figure that compares a
# thermistor, and neither needs to know which.
#
# Where a figure has a counterpart in METreport 11/2017 the docstring names it, so a FjordSim figure
# can be read beside the one that published model evaluation already drew.
#
# `src/Plotting.jl` is the pre-run counterpart of this file — it draws prepared *input* files. The
# shared idioms (`default_figure_size`, ivory land, finite-extrema colour ranges) come from there.

const VALIDATION_FIGURE_SIZE = (1100, 700)
const OBSERVED_COLOUR = :firebrick
const MODELLED_COLOUR = :black

"""
    finite_range(values...; pad = 0.05)

A colour or axis range over whatever is finite in `values`, padded, and `(0, 1)` when there is
nothing or everything is equal. The same guard `src/Plotting.jl` uses, which exists because a
station inside land is all `NaN` and Makie raises on an empty range rather than drawing nothing.
"""
function finite_range(values...; pad = 0.05)
    finite = Float64[]
    for collection in values
        append!(finite, (Float64(v) for v in collection if isfinite(v)))
    end

    isempty(finite) && return (0.0, 1.0)
    low, high = extrema(finite)
    low == high && return (low - 1, high + 1)
    margin = pad * (high - low)
    return (low - margin, high + margin)
end

"""
    plot_timeseries_comparison(path, observed, modelled; title, ylabel, constituents)

Observed against modelled at one station, in two panels: the tidal reconstruction of both, and the
residual each leaves behind.

METreport 11/2017 Fig. 8. The split is the whole point of that figure — the tide is the part a model
forced at an open boundary should get right almost for free, and the residual is the part that
actually tests the surge and the density field. Reporting only the total hides which of the two is
wrong.

Both series are fitted with the same constituents and the same phase reference, which is what makes
their amplitudes and phases comparable at all.
"""
function plot_timeseries_comparison(
    path,
    observed::ObservationSeries,
    modelled::ModelSeries;
    title = observed.station,
    ylabel = "$(observed.variable) ($(observed.units))",
    constituents = TIDAL_CONSTITUENTS,
    depth = nothing,
)
    times, o, m = match_series(observed, modelled; depth)
    isempty(times) && (@warn "No overlap to plot at $(observed.station)"; return nothing)

    reference = first(times)
    usable = [c for c in constituents if c.period_hours <= elapsed_hours(times)]
    if length(times) <= 1 + 2 * length(usable)
        @warn "Too short an overlap at $(observed.station) to fit a tide; skipping the figure." samples =
            length(times)
        return nothing
    end
    observed_fit = harmonic_fit(times, o; constituents = usable, reference)
    modelled_fit = harmonic_fit(times, m; constituents = usable, reference)

    hours = elapsed_hours.(times, reference)
    figure = Figure(size = VALIDATION_FIGURE_SIZE)

    tide_axis = Axis(
        figure[1, 1];
        title = "$title — combined tidal elevation",
        xlabel = "hours from $(reference)",
        ylabel,
    )
    lines!(tide_axis, hours, reconstruct(observed_fit, times); color = OBSERVED_COLOUR, label = "observed")
    lines!(tide_axis, hours, reconstruct(modelled_fit, times); color = MODELLED_COLOUR, label = "modelled")
    axislegend(tide_axis; position = :rt)

    residual_axis = Axis(
        figure[2, 1];
        title = "residual (total minus the fitted constituents)",
        xlabel = "hours from $(reference)",
        ylabel,
    )
    lines!(residual_axis, hours, o .- reconstruct(observed_fit, times); color = OBSERVED_COLOUR)
    lines!(residual_axis, hours, m .- reconstruct(modelled_fit, times); color = MODELLED_COLOUR)

    save(path, figure)
    return path
end

"""
    elapsed_hours(time, reference)
    elapsed_hours(times)

Hours from `reference` to `time`, or the span of `times`. The unit `TidalConstituent` periods are in.
"""
elapsed_hours(time::DateTime, reference::DateTime) =
    Dates.value(Dates.Millisecond(time - reference)) / 3_600_000
elapsed_hours(times) = isempty(times) ? 0.0 : elapsed_hours(last(times), first(times))

"""
    plot_tidal_constituents(path, fits; title)

Amplitude and phase per constituent, observed against modelled, with the standard errors the fits
report as whiskers.

METreport 11/2017 Table 3 as a figure. A table is the right medium for a few stations and thirteen
constituents, and `tidal_constituent_table` writes one; this is for seeing at a glance which
constituents the model has and which it invents — the report's own finding was that M2 came out well
everywhere while the shallow-water constituents M4, MN4 and MS4 came out with the right amplitude
and the wrong phase.

`fits` is `label => HarmonicFit`, with "observed" expected first.
"""
function plot_tidal_constituents(path, fits::AbstractVector{<:Pair}; title = "Tidal constituents")
    isempty(fits) && return nothing
    constituents = first(fits).second.constituents
    positions = 1:length(constituents)
    colours = [OBSERVED_COLOUR, MODELLED_COLOUR, :steelblue, :darkorange]

    figure = Figure(size = VALIDATION_FIGURE_SIZE)
    amplitude_axis = Axis(
        figure[1, 1];
        title = "$title — amplitude",
        ylabel = "amplitude",
        xticks = (positions, [c.name for c in constituents]),
    )
    phase_axis = Axis(
        figure[2, 1];
        title = "phase, relative to the shared reference epoch",
        ylabel = "phase (degrees)",
        xticks = (positions, [c.name for c in constituents]),
    )

    for (index, (label, fit)) in enumerate(fits)
        colour = colours[mod1(index, length(colours))]
        offset = 0.16 * (index - (1 + length(fits)) / 2)

        errorbars!(amplitude_axis, positions .+ offset, fit.amplitude, fit.amplitude_error; color = colour)
        scatter!(amplitude_axis, positions .+ offset, fit.amplitude; color = colour, label = String(label))
        errorbars!(phase_axis, positions .+ offset, fit.phase, fit.phase_error; color = colour)
        scatter!(phase_axis, positions .+ offset, fit.phase; color = colour)
    end

    axislegend(amplitude_axis; position = :rt)
    save(path, figure)
    return path
end

"""
    plot_hovmoller(path, series...; title, colormap, colorrange)

Time against depth, one panel per series, on a shared colour scale.

METreport 11/2017 Figs. 13-14 for currents and Figs. 25-28 for hydrography, all of which stack the
observed panel above the modelled one. The shared scale is what makes the stack readable and is the
reason this takes the series together rather than being called once each.

Depth runs downward, as every profile figure in both reports does.
"""
function plot_hovmoller(
    path,
    series...;
    title = "",
    colormap = :balance,
    colorrange = nothing,
    labels = nothing,
)
    isempty(series) && return nothing
    range = isnothing(colorrange) ? finite_range((s.values for s in series)...; pad = 0.0) : colorrange

    figure = Figure(size = (VALIDATION_FIGURE_SIZE[1], 260 * length(series) + 120))
    for (index, s) in enumerate(series)
        label = isnothing(labels) ? "$(s.station) $(s.variable)" : labels[index]
        axis = Axis(
            figure[index, 1];
            title = index == 1 ? "$title — $label" : label,
            xlabel = index == length(series) ? "time" : "",
            ylabel = "depth (m)",
            yreversed = true,
        )
        heatmap!(
            axis,
            elapsed_hours.(s.times, first(s.times)) ./ 24,
            s.depths,
            s.values;
            colormap,
            colorrange = range,
            nan_color = :ivory,
        )
        index == length(series) && (axis.xlabel = "days from $(first(s.times))")
    end

    Colorbar(figure[1:length(series), 2]; colormap, colorrange = range)
    save(path, figure)
    return path
end

"""
    plot_profiles(path, pairs; title, xlabel)

Observed and modelled profiles at the dates casts were taken — observed solid, modelled dashed, one
colour per date.

METreport 11/2017 Figs. 22-24, including the line style, which is theirs. This is the figure that
showed FjordOs' upper layer too thin in the outer fjord and its Drammensfjord deep water far too
fresh, and it is the most direct test of the stratification there is.

`pairs` is `date => (observed_values, modelled_values, depths)`.
"""
function plot_profiles(path, pairs; title = "", xlabel = "")
    isempty(pairs) && return nothing
    figure = Figure(size = (620, 760))
    axis = Axis(figure[1, 1]; title, xlabel, ylabel = "depth (m)", yreversed = true)
    colours = Makie.to_colormap(:viridis)

    for (index, (date, (observed, modelled, depths))) in enumerate(pairs)
        # Spread the casts across the colormap, so the legend reads as a season rather than as an
        # arbitrary set of colours.
        position = (index - 1) / max(1, length(pairs) - 1)
        colour = colours[clamp(round(Int, 1 + position * (length(colours) - 1)), 1, length(colours))]
        lines!(axis, observed, depths; color = colour, label = string(Date(date)))
        lines!(axis, modelled, depths; color = colour, linestyle = :dash)
    end

    axislegend(axis; position = :rb, nbanks = 2)
    save(path, figure)
    return path
end

"""
    plot_cast_profiles(path, pairs; title, spin_up_end = nothing, columns = 7)

Every cast at one station beside the model's profile nearest it in time, one panel per cast and
variable: a row per variable of `pairs` for each block of `columns` casts. Every panel spans the
same depths, the deepest the casts or the model column reach. A cast before `spin_up_end` is on a
grey background.

`pairs` is the station's `(observed, modelled, tolerance)` per variable. A cast with no model record
within `tolerance` is drawn alone.
"""
function plot_cast_profiles(path, pairs; title = "", spin_up_end = nothing, columns = 7)
    casts = sort(unique(time for (observed, _, _) in pairs for time in observed.times))
    isempty(casts) && return nothing

    deepest = maximum(max(deepest_depth(observed), deepest_depth(modelled)) for (observed, modelled, _) in pairs)
    per_row = min(length(casts), columns)
    blocks = collect(Iterators.partition(casts, per_row))
    rows = length(pairs)
    figure = Figure(size = (max(230 * per_row + 120, 900), 380 * rows * length(blocks) + 60))

    for (block_index, block) in enumerate(blocks), (column, time) in enumerate(block)
        for (row, (observed, modelled, tolerance)) in enumerate(pairs)
            axis = Axis(
                figure[rows * (block_index - 1) + row, column];
                title = row == 1 ? Dates.format(time, "yyyy-mm-dd") : "",
                xlabel = "$(observed.variable) ($(observed.units))",
                ylabel = column == 1 ? "depth (m)" : "",
                yreversed = true,
                xticks = LinearTicks(4),
            )
            record = nearest_record(modelled.times, time, tolerance)
            isnothing(record) || lines!(
                axis, modelled.values[record, :], modelled.depths;
                color = MODELLED_COLOUR, linewidth = 2.5, label = "modelled",
            )
            cast = findfirst(==(time), observed.times)
            isnothing(cast) || lines!(
                axis, observed.values[cast, :], observed.depths;
                color = OBSERVED_COLOUR, linewidth = 2, label = "observed",
            )
            # Deepest first: `ylims!` with the low end first sets `yreversed` back to false.
            ylims!(axis, deepest, 0)
            shade_spin_up!(axis, time, spin_up_end)
            block_index == 1 && row == 1 && column == 1 && axislegend(axis; position = :lb, labelsize = 10)
        end
    end

    Label(figure[0, :], title; fontsize = 18, font = :bold)
    save(path, figure)
    return path
end

"""
    deepest_depth(series)

The deepest level of an observed or modelled series holding a finite value, `0.0` for none.
"""
deepest_depth(series) =
    maximum((series.depths[k] for k in axes(series.values, 2) if any(isfinite, view(series.values, :, k))); init = 0.0)

shade_spin_up!(axis, time, ::Nothing) = nothing
shade_spin_up!(axis, time, spin_up_end::DateTime) = time < spin_up_end && (axis.backgroundcolor = :gray95)

# Marker shapes for the observation sources on the station map, in the order the setup names them.
const MAP_MARKERS = (:circle, :diamond, :rect, :pentagon, :hexagon, :star5)
# How many of the largest river outlets the station map names.
const MAP_RIVER_LABELS = 10
# Station names per line of a label on the station map.
const MAP_NAMES_PER_LINE = 4

"""
    plot_station_map(path, longitude, latitude, depth, sources, rivers, edges)

The stations of every observation source on the model depth, the river outlets, and the open
lateral `edges` as red lines along the domain's border.

`sources` is `(label, stations, scored)` per source, one marker shape each: a station whose name is
in `scored` is filled, one the run was not scored against is open. A source's stations in one model
cell are one marker, labelled with all their names (`group_stations`). `rivers` is
`(name, longitude, latitude, discharge)` per outlet; the outlets in one cell are one triangle, sized
by their summed mean discharge (`group_rivers`), the `MAP_RIVER_LABELS` largest named. `depth` is
positive down, `NaN` on land.
"""
function plot_station_map(path, longitude, latitude, depth, sources, rivers, edges)
    ratio = (last(longitude) - first(longitude)) * cosd(mean(latitude)) / (last(latitude) - first(latitude))
    figure = Figure(size = (1300, clamp(round(Int, 1100 / ratio) + 160, 600, 2400)), fontsize = 18)
    axis = Axis(figure[1, 1]; xlabel = "longitude (°E)", ylabel = "latitude (°N)",
        title = "Stations and river outlets on the model depth", aspect = AxisAspect(ratio))
    depth_map = heatmap!(axis, longitude, latitude, depth; colormap = :deep,
        colorrange = finite_range(depth; pad = 0.0), nan_color = RGBf(0.88, 0.87, 0.84))
    Colorbar(figure[1, 2], depth_map; label = "model depth (m)", height = Relative(0.5))

    # The domain's border, half a cell beyond the outermost cell centres.
    half_width, half_height = (longitude[2] - longitude[1]) / 2, (latitude[2] - latitude[1]) / 2
    border = (first(longitude) - half_width, last(longitude) + half_width,
        first(latitude) - half_height, last(latitude) + half_height)
    for (index, edge) in enumerate(edges)
        lines!(axis, edge_line(Val(edge), border)...; color = :red, linewidth = 7,
            label = index == 1 ? "open boundary" : nothing)
    end

    outlets = group_rivers(rivers)
    if !isempty(outlets)
        size_of(discharge) = 10 + 25 * sqrt(max(discharge, 0) / maximum(outlet[4] for outlet in outlets))
        scatter!(axis, [outlet[2] for outlet in outlets], [outlet[3] for outlet in outlets];
            marker = :utriangle, color = :dodgerblue, strokecolor = :black, strokewidth = 1,
            markersize = [size_of(outlet[4]) for outlet in outlets], label = "river outlet")
        for (name, x, y, discharge) in first(sort(outlets; by = outlet -> -outlet[4]), MAP_RIVER_LABELS)
            text!(axis, x, y; text = @sprintf(" %s (%.0f m³/s)", name, discharge), fontsize = 13,
                color = :navy, align = (:left, :bottom))
        end
    end

    for (index, (label, stations, scored)) in enumerate(sources)
        marker = MAP_MARKERS[mod1(index, length(MAP_MARKERS))]
        colour = Makie.wong_colors()[mod1(index, 7)]
        for (filled, suffix) in ((true, ""), (false, ", not scored"))
            groups = group_stations([station for station in stations if (station.name in scored) == filled], longitude, latitude)
            isempty(groups) && continue
            scatter!(axis, [group[2] for group in groups], [group[3] for group in groups];
                marker, markersize = 16, color = filled ? colour : :white, strokecolor = filled ? :black : colour,
                strokewidth = 2, label = label * suffix)
            for (names, x, y) in groups
                text!(axis, x, y; text = names, fontsize = 12, align = (:left, :center), font = :bold,
                    offset = (8, 0))
            end
        end
    end

    limits!(axis, border...)
    (isempty(rivers) && isempty(edges) && all(source -> isempty(source[2]), sources)) || axislegend(axis; position = :lt, framevisible = true, labelsize = 14)
    save(path, figure)
    return path
end

"""
    group_rivers(rivers)

The outlets of `rivers` merged by position, as `(names, longitude, latitude, discharge)`: the names
joined by `+` from the largest down and the discharges summed. `add_rivers` moves each outlet to the
head of the river bed the grid resolves, so the rivers of one fjord head share one position, though
each keeps inflow columns of its own.
"""
function group_rivers(rivers)
    positions = unique((river[2], river[3]) for river in rivers)
    return map(positions) do (x, y)
        members = sort([river for river in rivers if (river[2], river[3]) == (x, y)]; by = river -> -river[4])
        (join(first.(members), " + "), x, y, sum(river[4] for river in members))
    end
end

"""
    group_stations(stations, longitude, latitude)

`stations` merged by the model cell nearest each, as `(names, longitude, latitude)`: the names in
the order given, `MAP_NAMES_PER_LINE` to a line, at the stations' mean position.
"""
function group_stations(stations, longitude, latitude)
    cell(station) = (argmin(abs.(longitude .- station.longitude)), argmin(abs.(latitude .- station.latitude)))
    cells = unique(cell.(stations))
    return map(cells) do position
        members = [station for station in stations if cell(station) == position]
        names = join([join(line, ", ") for line in Iterators.partition([station.name for station in members], MAP_NAMES_PER_LINE)], "\n")
        (names, mean(station.longitude for station in members), mean(station.latitude for station in members))
    end
end

"""
    edge_line(edge, (west, east, south, north))

The two ends of one lateral edge of the border, as `(longitudes, latitudes)`.
"""
edge_line(::Val{:south}, (west, east, south, north)) = ([west, east], [south, south])
edge_line(::Val{:north}, (west, east, south, north)) = ([west, east], [north, north])
edge_line(::Val{:west}, (west, east, south, north)) = ([west, west], [south, north])
edge_line(::Val{:east}, (west, east, south, north)) = ([east, east], [south, north])

"""
    plot_drift(path, entries, start, stop; variables = ("S", "T"))

Each station at one depth through `start` to `stop`: the model at every record as a line, the casts
on their days as points. One row per station, one column per variable.

`entries` is `(station, depth, pairs)`, `pairs` the station's `(observed, modelled, tolerance)` per
variable. A cast's value is the mean of its samples within 1 m of `depth`; the model's is its
profile interpolated to `depth`.
"""
function plot_drift(path, entries, start, stop; variables = ("S", "T"))
    isempty(entries) && return nothing
    days(time) = elapsed_hours(time, start) / 24

    figure = Figure(size = (1500, 270 * length(entries) + 40))
    for (row, (station, depth, pairs)) in enumerate(entries), (column, variable) in enumerate(variables)
        index = findfirst(pair -> pair[1].variable == variable, pairs)
        isnothing(index) && continue
        observed, modelled, _ = pairs[index]

        axis = Axis(
            figure[row, column];
            title = @sprintf("%s, %s at %.0f m", station, variable, depth),
            ylabel = observed.units,
            xlabel = row == length(entries) ? "days since $(Date(start))" : "",
        )
        modelled_values = [first(profile_at(modelled.depths, view(modelled.values, n, :), [depth])) for n in eachindex(modelled.times)]
        lines!(axis, days.(modelled.times), modelled_values; color = MODELLED_COLOUR, linewidth = 2, label = "modelled")

        near = [k for k in eachindex(observed.depths) if abs(observed.depths[k] - depth) <= 1]
        observed_values = [cast_mean(view(observed.values, n, near)) for n in eachindex(observed.times)]
        scatter!(axis, days.(observed.times), observed_values; color = OBSERVED_COLOUR, markersize = 12, label = "observed")

        xlims!(axis, 0, days(stop))
        row == 1 && column == 1 && axislegend(axis; position = :lb, labelsize = 10)
    end

    save(path, figure)
    return path
end

cast_mean(values) = (finite = filter(isfinite, values); isempty(finite) ? NaN : mean(finite))

"""
    plot_current_rose(path, observed, modelled; title, bins)

Speed and direction distributions, observed beside modelled: a direction histogram on a polar axis
and a speed histogram beneath it.

METreport 11/2017 Figs. 19 and 20. `observed` and `modelled` are each `(eastward, northward)`.
"""
function plot_current_rose(path, observed, modelled; title = "", bins = 36)
    figure = Figure(size = VALIDATION_FIGURE_SIZE)

    for (column, (label, components, colour)) in enumerate((
        ("observed", observed, OBSERVED_COLOUR),
        ("modelled", modelled, MODELLED_COLOUR),
    ))
        speed, direction = direction_statistics(components[1], components[2])
        finite = findall(i -> isfinite(speed[i]) && isfinite(direction[i]), eachindex(speed))

        rose = PolarAxis(figure[1, column]; title = "$title — $label", theta_0 = π / 2, direction = -1)
        if !isempty(finite)
            edges = range(0, 2π; length = bins + 1)
            counts = zeros(Int, bins)
            for index in finite
                counts[clamp(floor(Int, direction[index] / 360 * bins) + 1, 1, bins)] += 1
            end
            centres = [(edges[k] + edges[k + 1]) / 2 for k = 1:bins]
            barplot!(rose, centres, counts ./ sum(counts); width = 2π / bins, color = colour)
        end

        histogram = Axis(figure[2, column]; xlabel = "speed (m/s)", ylabel = "density")
        isempty(finite) || hist!(histogram, speed[finite]; bins = 60, normalization = :pdf, color = colour)
    end

    save(path, figure)
    return path
end

"""
    plot_taylor(path, entries; title)

A Taylor diagram: every station-variable pair as one point, placed by its correlation (as the polar
angle) and its normalised standard deviation (as the radius), so that its distance from the
reference point on the axis is its normalised centred RMSE.

The one figure that answers "how does this run score overall" at a glance, and a summary the FjordOs
reports predate. A point on the arc at radius 1 has the right variability, a point near the axis has
the right phasing, and the reference marker at `(1, 0)` is a perfect model.

`entries` is `label => SkillMetrics`.
"""
function plot_taylor(path, entries::AbstractVector{<:Pair}; title = "Taylor diagram")
    figure = Figure(size = (760, 760))
    axis = PolarAxis(
        figure[1, 1];
        title,
        thetalimits = (0, π / 2),
        rticks = 0:0.5:2,
        thetaticks = (acos.([1.0, 0.99, 0.95, 0.9, 0.8, 0.6, 0.4, 0.2, 0.0]),
                      ["1", "0.99", "0.95", "0.9", "0.8", "0.6", "0.4", "0.2", "0"]),
    )

    # The reference: perfect correlation, observed variability.
    scatter!(axis, [0.0], [1.0]; marker = :star5, markersize = 22, color = :grey25)

    colours = Makie.to_colormap(:tab20)
    elements = Makie.MarkerElement[]
    labels = String[]

    for (index, (label, metrics)) in enumerate(entries)
        isfinite(metrics.correlation) && isfinite(metrics.std_ratio) || continue
        colour = colours[mod1(index, length(colours))]
        scatter!(
            axis,
            [acos(clamp(metrics.correlation, -1, 1))],
            [metrics.std_ratio];
            color = colour,
            markersize = 13,
        )
        push!(elements, MarkerElement(color = colour, marker = :circle, markersize = 13))
        push!(labels, String(label))
    end

    # Built from explicit elements rather than by handing `Legend` the axis: Makie's
    # `get_labeled_plots` cannot scrape a `PolarAxis`, so the usual `Legend(figure[1, 2], axis)`
    # raises here where it works for a Cartesian one.
    isempty(elements) || Legend(figure[1, 2], elements, labels; framevisible = false)
    save(path, figure)
    return path
end

"""
    plot_target(path, entries; title)

A target diagram: normalised bias against normalised, sign-carrying centred RMSE, so that a point's
distance from the origin is its normalised total RMSE.

The complement to `plot_taylor`, which cannot show bias at all. Inside the unit circle the model
beats predicting the observed mean; left of the axis it is under-energetic, right of it
over-energetic — which is the axis the reports' repeated "simulated currents are stronger than the
observed" belongs on.
"""
function plot_target(path, entries::AbstractVector{<:Pair}; title = "Target diagram")
    figure = Figure(size = (760, 720))
    axis = Axis(
        figure[1, 1];
        title,
        xlabel = "sign(σ_model − σ_obs) × centred RMSE / σ_obs",
        ylabel = "bias / σ_obs",
        aspect = 1,
    )

    angles = range(0, 2π; length = 200)
    lines!(axis, cos.(angles), sin.(angles); color = :grey, linestyle = :dash)
    hlines!(axis, [0.0]; color = :grey85)
    vlines!(axis, [0.0]; color = :grey85)

    colours = Makie.to_colormap(:tab20)
    drawn = false
    for (index, (label, metrics)) in enumerate(entries)
        isfinite(metrics.centred_rmse) && isfinite(metrics.bias) && metrics.observed_std > 0 || continue
        drawn = true
        sign_of = metrics.modelled_std >= metrics.observed_std ? 1 : -1
        scatter!(
            axis,
            [sign_of * metrics.centred_rmse / metrics.observed_std],
            [metrics.bias / metrics.observed_std];
            color = colours[mod1(index, length(colours))],
            markersize = 13,
            label = String(label),
        )
    end

    # `Legend` raises on an axis with no labelled plot.
    drawn && Legend(figure[1, 2], axis; framevisible = false)
    save(path, figure)
    return path
end

"""
    skill_table(entries)

The scores as a plain-text table, one row per station-variable pair.

A figure is the wrong medium for numbers a reader wants to quote, and both reports give their skill
summaries as tables for that reason.
"""
function skill_table(entries::AbstractVector{<:Pair})
    header = @sprintf(
        "%-28s %7s %9s %9s %9s %9s %7s %7s %7s\n",
        "station / variable", "n", "obs mean", "bias", "RMSE", "cRMSE", "r", "σm/σo", "d",
    )
    rule = repeat("-", length(header) - 1) * "\n"
    rows = map(entries) do (label, m)
        @sprintf(
            "%-28s %7d %9.3f %9.3f %9.3f %9.3f %7.3f %7.3f %7.3f\n",
            String(label), m.count, m.observed_mean, m.bias, m.rmse,
            m.centred_rmse, m.correlation, m.std_ratio, m.willmott,
        )
    end

    return string(header, rule, rows...)
end

"""
    tidal_constituent_table(fits)

Amplitude and phase per constituent for each fit, as a plain-text table in the shape of METreport
11/2017 Table 3 — constituent and period down the side, each series' amplitude and phase across.
"""
function tidal_constituent_table(fits::AbstractVector{<:Pair})
    isempty(fits) && return ""
    constituents = first(fits).second.constituents

    header = @sprintf("%-6s %9s", "comp.", "period[h]")
    for (label, _) in fits
        header *= @sprintf(" %10s %8s", "$(label) amp", "phase")
    end
    header *= "\n"

    rows = map(enumerate(constituents)) do (index, constituent)
        row = @sprintf("%-6s %9.4f", constituent.name, constituent.period_hours)
        for (_, fit) in fits
            row *= @sprintf(" %10.4f %8.1f", fit.amplitude[index], fit.phase[index])
        end
        row * "\n"
    end

    return string(header, repeat("-", length(header) - 1), "\n", rows...)
end

"""
    plot_section(path, snapshot_path, variable; latitude, record = 1, colormap = :balance, title)

A vertical section across the fjord at constant latitude, with the model's own bathymetry drawn over
it.

METreport 11/2017 Fig. 11 draws its current transects this way. A transect that is very nearly
zonal is one row of a lat-lon grid and needs no interpolation, which is why this takes a latitude
rather than two endpoints.

Read from a `SnapshotWriter` file rather than from the station writers, because a section is the one
diagnostic that wants the whole field and does not want sub-daily sampling. The bathymetry comes
from the same file — `bottom_height` is written beside the fields — so the section and the seabed
cannot disagree.

That report's version also overlays the *real* bathymetry to show how far the model's had to be
smoothed, which is its main point there. That comparison is not drawn here because it does not
apply: FjordSim's z-coordinate needs no rx0 smoothing.
"""
function plot_section(
    path,
    snapshot_path,
    variable;
    latitude,
    record = 1,
    colormap = :balance,
    title = "",
)
    return NCDataset(snapshot_path) do ds
        haskey(ds, variable) || throw(
            ArgumentError(
                "Snapshot $snapshot_path has no variable $variable. It carries " *
                "$(join(sort(collect(keys(ds))), ", ")).",
            ),
        )

        raw = ds[variable]
        dimensions = collect(NCDatasets.dimnames(raw))
        # Named by field location, so a Face-located velocity and a Center-located tracer do not
        # share dimension names; pick them out by prefix rather than by position.
        longitude_name = dimensions[findfirst(name -> startswith(name, "λ"), dimensions)]
        latitude_name = dimensions[findfirst(name -> startswith(name, "φ"), dimensions)]
        vertical_name = dimensions[findfirst(name -> startswith(name, "z"), dimensions)]

        longitudes = Array{Float64}(ds[longitude_name][:])
        latitudes = Array{Float64}(ds[latitude_name][:])
        depths = Array{Float64}(ds[vertical_name][:])

        row = argmin(abs.(latitudes .- latitude))
        slab = Array(raw[:, row, :, record])
        values = [ismissing(v) ? NaN : Float64(v) for v in slab]

        figure = Figure(size = VALIDATION_FIGURE_SIZE)
        axis = Axis(
            figure[1, 1];
            title = "$title — $variable at $(round(latitudes[row]; digits = 4))°N",
            xlabel = "longitude (°E)",
            ylabel = "depth (m)",
        )

        range = finite_range(values; pad = 0.0)
        heatmap!(axis, longitudes, depths, values; colormap, colorrange = range, nan_color = :ivory)
        # Zero, which on a current section is the line between inflow and outflow — the contour
        # METreport 11/2017 Fig. 11 draws in black.
        any(v -> v < 0, filter(isfinite, values)) &&
            any(v -> v > 0, filter(isfinite, values)) &&
            contour!(axis, longitudes, depths, values; levels = [0.0], color = :black)

        if haskey(ds, "bottom_height")
            bottom = Array{Float64}(ds["bottom_height"][:, row, 1])
            # `bottom_height` is on tracer columns; a Face-located variable has one more of them.
            lines!(axis, longitudes[1:min(end, length(bottom))], bottom[1:min(end, length(longitudes))];
                   color = :black, linewidth = 2.5)
        end

        Colorbar(figure[1, 2]; colormap, colorrange = range)
        save(path, figure)
        return path
    end
end

"""
    plot_domain_extremes(path, statistics, ranges)

The domain's minimum, volume mean and maximum of each variable through the run, one row each, with
the plausible range `ranges` states dashed; and beside them how many wet cells lie outside it. A
runaway cell shows as a spike in the extremes long before it moves the mean.
"""
function plot_domain_extremes(path, statistics, ranges)
    variables = sort(unique(s.variable for s in statistics))
    figure = Figure(size = (VALIDATION_FIGURE_SIZE[1], 320 * length(variables)))

    for (row, variable) in enumerate(variables)
        rows = [s for s in statistics if s.variable == variable]
        days = [Dates.value(s.time - first(rows).time) / 86_400_000 for s in rows]
        xlabel = "days since $(Date(first(rows).time))"

        axis = Axis(figure[row, 1]; title = "$variable over the domain", xlabel)
        lines!(axis, days, [s.minimum for s in rows]; color = :dodgerblue, label = "minimum")
        lines!(axis, days, [s.mean for s in rows]; color = MODELLED_COLOUR, label = "volume mean")
        lines!(axis, days, [s.maximum for s in rows]; color = OBSERVED_COLOUR, label = "maximum")
        if haskey(ranges, variable)
            hlines!(axis, collect(ranges[variable]); color = :grey, linestyle = :dash)
        end
        axislegend(axis; position = :rb, labelsize = 11)

        counts = Axis(figure[row, 2]; title = "wet cells outside the range", xlabel)
        lines!(counts, days, [s.below for s in rows]; color = :dodgerblue, label = "below")
        lines!(counts, days, [s.above for s in rows]; color = OBSERVED_COLOUR, label = "above")
        axislegend(counts; position = :rt, labelsize = 11)
    end

    save(path, figure)
    return path
end
