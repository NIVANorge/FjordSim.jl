# Validation figures. Every function takes series and metrics rather than configs or paths, for the
# same reason `Metrics.jl` does: the figure that compares a tide gauge is the figure that compares a
# thermistor, and neither needs to know which.
#
# Where a figure has a counterpart in METreport 11/2017 the docstring names it, because the point of
# the exercise is to put FjordSim's answer beside the one that report already published.
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

METreport 11/2017 Table 3 as a figure. A table is the right medium for three stations and thirteen
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
    plot_qq_scatter(path, observed, modelled; title, unit)

A quantile-quantile plot beside a scatter plot of the same pairs.

METreport 11/2017 Figs. 21 and 30, and the pair is deliberate. The QQ plot asks whether the model
has the right *distribution* and the scatter plot whether it has it at the right *time*; the report's
own conclusion at Slagen was that the first was nearly right and the second was not, which is a very
different verdict from either panel alone.
"""
function plot_qq_scatter(path, observed, modelled; title = "", unit = "")
    probabilities = 0.005:0.005:0.995
    observed_quantiles = quantiles(observed, probabilities)
    modelled_quantiles = quantiles(modelled, probabilities)
    limits = finite_range(observed, modelled)

    figure = Figure(size = VALIDATION_FIGURE_SIZE)
    label = isempty(unit) ? "" : " ($unit)"

    qq = Axis(
        figure[1, 1];
        title = "$title — quantile-quantile",
        xlabel = "observed$label",
        ylabel = "modelled$label",
        aspect = 1,
    )
    lines!(qq, [limits...], [limits...]; color = :grey, linestyle = :dash)
    scatter!(qq, observed_quantiles, modelled_quantiles; color = MODELLED_COLOUR, markersize = 5)
    limits!(qq, limits..., limits...)

    scatter_axis = Axis(
        figure[1, 2];
        title = "pairs in time",
        xlabel = "observed$label",
        ylabel = "modelled$label",
        aspect = 1,
    )
    lines!(scatter_axis, [limits...], [limits...]; color = :grey, linestyle = :dash)
    scatter!(scatter_axis, observed, modelled; color = (:steelblue, 0.25), markersize = 4)
    limits!(scatter_axis, limits..., limits...)

    save(path, figure)
    return path
end

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

The one figure that answers "is this comparable to FjordOs" at a glance, and the summary both
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
    for (index, (label, metrics)) in enumerate(entries)
        isfinite(metrics.centred_rmse) && isfinite(metrics.bias) && metrics.observed_std > 0 || continue
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

    Legend(figure[1, 2], axis; framevisible = false)
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

METreport 11/2017 Fig. 11, which is how both Statnett transects are shown: the Filtvedt-Brenntangen
line at 59.582°N and the Småskjær-Evje line at ~59.35°N are both very nearly zonal, so each is one
row of a lat-lon grid and needs no interpolation.

Read from a `SnapshotWriter` file rather than from the station writers, because a section is the one
diagnostic that wants the whole field and does not want sub-daily sampling. The bathymetry comes
from the same file — `bottom_height` is written beside the fields — so the section and the seabed
cannot disagree.

That report's version also overlays the *real* bathymetry to show how far the model's had to be
smoothed, which is its main point there. That comparison is not drawn here because it does not
apply: FjordSim's z-coordinate needs no rx0 smoothing, and `max_slope_factor = 0.25` costs 0.6 m at
the Drøbak sill against a 395.1 m basin.
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
