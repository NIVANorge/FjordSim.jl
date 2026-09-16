# Model-observation statistics. Nothing here knows what a fjord is or where a file lives: every
# function takes paired vectors, or a time axis and a series, and returns numbers. That is what lets
# the same code score a tide gauge, a current meter and a CTD cast.
#
# The set is the one the two FjordOs reports use, plus the summary diagrams they predate. METreport
# 11/2017 scores with monthly means and variances (Table 6), harmonic amplitude and phase (Table 3),
# a 49-hour running mean to separate the estuarine circulation (its eq. 3), the barotropic/baroclinic
# split (its eqs. 1-2), and QQ and scatter plots; Taylor and target diagrams are the standard way to
# put all of that on one axis and are built from `SkillMetrics` below.

"""
    SkillMetrics

The summary of one model-observation comparison. Built by `skill_metrics`.

# Fields
- `count`: pairs that survived the finite filter.
- `observed_mean`, `modelled_mean`, `observed_std`, `modelled_std`: the two distributions.
- `bias`: `modelled_mean - observed_mean`. Positive means the model runs high.
- `rmse`: root mean square error, bias included.
- `centred_rmse`: the same with the bias removed, so `rmse² = bias² + centred_rmse²`. This is the
  radial coordinate of a Taylor diagram and the one to quote when a constant offset is not the
  interesting part of the error.
- `correlation`: Pearson correlation in time. The reports are blunt that this is the metric a
  high-resolution model does worst on, and that the reason is physical rather than a defect: a model
  resolving eddies and meanders puts them in the right statistics at the wrong moment.
- `std_ratio`: `modelled_std / observed_std`. Above one means the model is too energetic.
- `willmott`: Willmott's index of agreement, `1 - Σ(m-o)² / Σ(|m-ō| + |o-ō|)²`, in `[0, 1]` with 1
  perfect. Bounded and defined even when the observations barely vary, which is why it is preferred
  to `murphy` for a short record.
- `murphy`: Murphy skill score, `1 - MSE / Var(o)`, identical to Nash-Sutcliffe. Zero means the model
  is no better than predicting the observed mean, and negative means worse; unbounded below.
"""
struct SkillMetrics
    count::Int
    observed_mean::Float64
    modelled_mean::Float64
    observed_std::Float64
    modelled_std::Float64
    bias::Float64
    rmse::Float64
    centred_rmse::Float64
    correlation::Float64
    std_ratio::Float64
    willmott::Float64
    murphy::Float64
end

"""
    skill_metrics(observed, modelled)

Score two aligned series against each other. Pairs where either value is not finite are dropped, so
a masked model cell or a gap in an instrument record costs only those samples.

Returns a `SkillMetrics` of `NaN`s with `count = 0` when nothing survives, rather than raising: a
station with no overlap is a normal outcome of a validation sweep and should show up as a blank row,
not stop the run.
"""
function skill_metrics(observed, modelled)
    length(observed) == length(modelled) || throw(
        DimensionMismatch(
            "skill_metrics needs aligned series, got $(length(observed)) observations and " *
            "$(length(modelled)) model values.",
        ),
    )

    keep = [i for i in eachindex(observed) if isfinite(observed[i]) && isfinite(modelled[i])]
    length(keep) >= 2 || return SkillMetrics(length(keep), ntuple(_ -> NaN, 11)...)

    o = Float64[observed[i] for i in keep]
    m = Float64[modelled[i] for i in keep]

    ō = mean(o)
    m̄ = mean(m)
    σo = std(o; corrected = false)
    σm = std(m; corrected = false)

    bias = m̄ - ō
    mse = mean((m .- o) .^ 2)
    rmse = sqrt(mse)
    centred = sqrt(max(zero(mse), mse - bias^2))

    # `cor` is undefined for a constant series; a flat model or a flat instrument is a real outcome.
    correlation = σo > 0 && σm > 0 ? cor(o, m) : NaN

    denominator = sum((abs.(m .- ō) .+ abs.(o .- ō)) .^ 2)
    willmott = denominator > 0 ? 1 - sum((m .- o) .^ 2) / denominator : NaN
    murphy = σo > 0 ? 1 - mse / σo^2 : NaN

    return SkillMetrics(
        length(keep), ō, m̄, σo, σm, bias, rmse, centred, correlation,
        σo > 0 ? σm / σo : NaN, willmott, murphy,
    )
end

"""
    monthly_statistics(times, values)

Count, mean and variance of `values` per calendar month, as a `Dict` keyed by `(year, month)`.

METreport 11/2017 Table 6 reports exactly this for observed and simulated temperature at
Åsgårdstrand, and the count matters as much as the moments there: the observed and simulated counts
differ by a third in every row because the instrument has gaps and the model does not, so a mean
difference cannot be read without knowing how many samples each came from.
"""
function monthly_statistics(times, values)
    buckets = Dict{Tuple{Int,Int},Vector{Float64}}()
    for (time, value) in zip(times, values)
        isfinite(value) || continue
        push!(get!(buckets, (Dates.year(time), Dates.month(time)), Float64[]), value)
    end

    return Dict(
        key => (count = length(v), mean = mean(v), variance = var(v; corrected = false))
        for (key, v) in buckets
    )
end

"""
    running_mean(times, values, window)

Centred running mean of `values` over a `window` of time, as a vector the same length as `values`,
with `NaN` where the window is not fully covered.

`window` is a `Period` rather than a sample count so the same call works on a 10-minute tide-gauge
series and an hourly current record. METreport 11/2017 uses `Hour(49)`, which is the conventional
choice for a fjord: it is four M2 cycles and two inertial periods at this latitude, so it removes
the tide without removing the estuarine signal underneath it.
"""
function running_mean(times, values, window::Dates.Period)
    half = Dates.Millisecond(window) ÷ 2
    output = fill(NaN, length(values))
    isempty(values) && return output

    first_time, last_time = first(times), last(times)
    for i in eachindex(values)
        centre = times[i]
        centre - half >= first_time && centre + half <= last_time || continue

        total = 0.0
        count = 0
        for j in eachindex(values)
            times[j] < centre - half && continue
            times[j] > centre + half && break
            isfinite(values[j]) || continue
            total += values[j]
            count += 1
        end

        count > 0 && (output[i] = total / count)
    end

    return output
end

"""
    depth_average(profile, thicknesses)

`u₀ = (1/h) ∫ u dz`, the barotropic mode of METreport 11/2017 eq. 1, over whatever part of the
column is wet.

`profile` and `thicknesses` run bottom to top as the model's own vertical index does. Dry cells are
those whose value is not finite — the same land sentinel `ForcingFromFile` uses — so a partial
bottom cell contributes its real thickness and rock contributes nothing.
"""
function depth_average(profile, thicknesses)
    total = 0.0
    depth = 0.0
    for k in eachindex(profile)
        isfinite(profile[k]) || continue
        total += profile[k] * thicknesses[k]
        depth += thicknesses[k]
    end

    return depth > 0 ? total / depth : NaN
end

"""
    baroclinic_deviation(profile, thicknesses)

`uₙ = u - u₀`, the baroclinic mode of METreport 11/2017 eq. 2: the profile with its depth average
removed. Dry cells stay `NaN`.
"""
function baroclinic_deviation(profile, thicknesses)
    mean_value = depth_average(profile, thicknesses)
    return [isfinite(value) ? value - mean_value : NaN for value in profile]
end

# --- Harmonic analysis ---

"""
    TidalConstituent

A tidal constituent, named and with its period in hours.
"""
struct TidalConstituent
    name::String
    period_hours::Float64
end

"""
    TIDAL_CONSTITUENTS

The thirteen constituents METreport 11/2017 Table 3 reports, in that table's order (by period,
longest first). Periods are the standard astronomical ones the report quotes.

Only eleven of them were in FjordOs' own tidal forcing; SA and SSA were not, and the report's point
in listing them anyway is that the model picks them up regardless, through the daily-mean sea level
it takes from NorKyst at the open boundary. That test applies here unchanged — and more directly,
since this setup imposes no harmonic forcing at all and takes its entire tide from the boundary
data.
"""
const TIDAL_CONSTITUENTS = [
    TidalConstituent("SA", 8764.0),
    TidalConstituent("SSA", 4382.0),
    TidalConstituent("Q1", 26.8684),
    TidalConstituent("O1", 25.8193),
    TidalConstituent("P1", 24.0659),
    TidalConstituent("K1", 23.9345),
    TidalConstituent("N2", 12.6584),
    TidalConstituent("M2", 12.4206),
    TidalConstituent("S2", 12.0000),
    TidalConstituent("K2", 11.9672),
    TidalConstituent("MN4", 6.2692),
    TidalConstituent("M4", 6.2103),
    TidalConstituent("MS4", 6.1033),
]

"""
    HarmonicFit

The result of `harmonic_fit`: a mean, and an amplitude and phase per constituent.

# Fields
- `constituents`: what was fitted, in the order the other vectors are indexed by.
- `mean`: the constant term, in the units of the input.
- `amplitude`, `phase`: per constituent; `phase` in degrees in `[0, 360)`, measured from `reference`.
- `amplitude_error`, `phase_error`: one standard error each, from the residual variance.
- `reference`: the epoch phases are measured from.
- `residual_std`: standard deviation of what the fit did not explain — the non-tidal residual, which
  is the lower panel of METreport 11/2017 Fig. 8.
- `count`: finite samples the fit used.

**Phases are relative to `reference`, not to Greenwich.** No nodal correction or astronomical
argument is applied, so a phase here is not directly comparable to a published tide-table phase.
What *is* comparable, and is what a validation needs, is the difference between two fits made with
the same `reference` — so score model against observations by fitting both and differencing, exactly
as METreport 11/2017 does with `t_tide`.
"""
struct HarmonicFit
    constituents::Vector{TidalConstituent}
    mean::Float64
    amplitude::Vector{Float64}
    phase::Vector{Float64}
    amplitude_error::Vector{Float64}
    phase_error::Vector{Float64}
    reference::DateTime
    residual_std::Float64
    count::Int
end

"""
    harmonic_fit(times, values; constituents = TIDAL_CONSTITUENTS, reference = first(times))

Least-squares harmonic analysis: fit `mean + Σ aₖcos(ωₖt) + bₖsin(ωₖt)` and return amplitudes and
phases.

This is what `t_tide` does for METreport 11/2017, minus the nodal corrections — see `HarmonicFit` on
why that does not affect a model-observation comparison. Non-finite samples are dropped, so a gappy
instrument record needs no pre-filling; unequal spacing is fine too, since nothing here assumes a
sampling interval.

A constituent is only resolvable if the record is long enough to separate it from its neighbours, and
a fit will happily return a number where it is not. Two constituents of periods `T₁` and `T₂` need a
record of about `T₁T₂/|T₁-T₂|` to separate — 183 days for S2 against K2, which is the tightest pair
in `TIDAL_CONSTITUENTS`, and about a year for SA against the mean. A warning is raised when the
record is shorter than the longest period fitted, which is the crudest form of that check and the
one that catches the mistake that actually happens.
"""
function harmonic_fit(
    times,
    values;
    constituents = TIDAL_CONSTITUENTS,
    reference = first(times),
)
    keep = [i for i in eachindex(values) if isfinite(values[i])]
    n = length(keep)
    m = 1 + 2 * length(constituents)

    n > m || throw(
        ArgumentError(
            "harmonic_fit needs more finite samples than unknowns: $n samples against $m unknowns " *
            "for $(length(constituents)) constituents.",
        ),
    )

    span_hours = Dates.value(Dates.Millisecond(times[keep[end]] - times[keep[1]])) / 3_600_000
    longest = maximum(constituent.period_hours for constituent in constituents)
    span_hours >= longest || @warn(
        "Record is shorter than the longest constituent fitted; its amplitude and phase are not " *
        "resolvable and the fit will absorb other signal into it.",
        record_days = round(span_hours / 24; digits = 1),
        longest_period_days = round(longest / 24; digits = 1),
    )

    # Hours since `reference`, which is what the periods are in.
    t = [Dates.value(Dates.Millisecond(times[i] - reference)) / 3_600_000 for i in keep]
    y = Float64[values[i] for i in keep]

    design = ones(Float64, n, m)
    for (k, constituent) in enumerate(constituents)
        ω = 2π / constituent.period_hours
        design[:, 2k] = cos.(ω .* t)
        design[:, 2k + 1] = sin.(ω .* t)
    end

    coefficients = design \ y
    residual = y .- design * coefficients
    # The unbiased residual variance, and the parameter covariance it implies.
    variance = sum(abs2, residual) / (n - m)
    covariance = variance .* inv(design' * design)

    amplitude = Float64[]
    phase = Float64[]
    amplitude_error = Float64[]
    phase_error = Float64[]

    for k in eachindex(constituents)
        a = coefficients[2k]
        b = coefficients[2k + 1]
        magnitude = hypot(a, b)
        push!(amplitude, magnitude)
        push!(phase, mod(rad2deg(atan(b, a)), 360))

        # Propagate the two coefficient variances through `hypot` and `atan`. With `a` and `b`
        # near-orthogonal in the design — which they are for any record covering several periods —
        # the cross term is small and the two errors reduce to the usual scaling.
        σa² = covariance[2k, 2k]
        σb² = covariance[2k + 1, 2k + 1]
        if magnitude > 0
            push!(amplitude_error, sqrt((a^2 * σa² + b^2 * σb²) / magnitude^2))
            push!(phase_error, rad2deg(sqrt((b^2 * σa² + a^2 * σb²) / magnitude^4)))
        else
            push!(amplitude_error, sqrt(max(σa², σb²)))
            push!(phase_error, NaN)
        end
    end

    return HarmonicFit(
        collect(constituents), coefficients[1], amplitude, phase,
        amplitude_error, phase_error, reference, std(residual; corrected = false), n,
    )
end

"""
    reconstruct(fit, times)

The elevation `fit` predicts at `times` — the "combined water elevation" of METreport 11/2017 §4.1,
which is the sum of the fitted constituents. Subtracting it from the original series gives that
section's "residual", the part the tide does not explain.
"""
function reconstruct(fit::HarmonicFit, times)
    output = fill(fit.mean, length(times))
    for (index, time) in enumerate(times)
        t = Dates.value(Dates.Millisecond(time - fit.reference)) / 3_600_000
        for k in eachindex(fit.constituents)
            ω = 2π / fit.constituents[k].period_hours
            output[index] += fit.amplitude[k] * cos(ω * t - deg2rad(fit.phase[k]))
        end
    end

    return output
end

"""
    TidalEllipse

One constituent's current ellipse: what a rotating tidal current traces out over a period.

# Fields
- `constituent`
- `semi_major`, `semi_minor`: axis lengths in the velocity's units. `semi_minor` is signed — positive
  anticlockwise, negative clockwise — which is the convention that keeps the sense of rotation.
- `inclination`: degrees anticlockwise from east to the major axis, in `[0, 180)`.
- `phase`: degrees from the fit's reference to maximum flow along the major axis.

METreport 11/2017 Table 4 reports the major amplitude and phase of the barotropic current at Km1
this way.
"""
struct TidalEllipse
    constituent::TidalConstituent
    semi_major::Float64
    semi_minor::Float64
    inclination::Float64
    phase::Float64
end

"""
    tidal_ellipses(eastward_fit, northward_fit)

Combine two `harmonic_fit`s of the same constituents — one per velocity component — into one ellipse
per constituent.

Via the rotary decomposition: a component pair at one frequency is the sum of an anticlockwise
circular current of radius `Wp` and a clockwise one of radius `Wm`, so the ellipse has semi-axes
`Wp ± Wm` and its orientation and phase are the half-sum and half-difference of the two rotary
phases. That is the standard construction and the one `t_tide` uses.
"""
function tidal_ellipses(eastward::HarmonicFit, northward::HarmonicFit)
    eastward.constituents == northward.constituents || throw(
        ArgumentError("tidal_ellipses needs two fits over the same constituents, in the same order."),
    )
    eastward.reference == northward.reference || throw(
        ArgumentError("tidal_ellipses needs two fits sharing one phase reference."),
    )

    ellipses = TidalEllipse[]
    for k in eachindex(eastward.constituents)
        θu = deg2rad(eastward.phase[k])
        θv = deg2rad(northward.phase[k])
        au = eastward.amplitude[k]
        av = northward.amplitude[k]

        # Rotary amplitudes and phases.
        wp = hypot(au * cos(θu) + av * sin(θv), au * sin(θu) - av * cos(θv)) / 2
        wm = hypot(au * cos(θu) - av * sin(θv), au * sin(θu) + av * cos(θv)) / 2
        θp = atan(au * sin(θu) - av * cos(θv), au * cos(θu) + av * sin(θv))
        θm = atan(au * sin(θu) + av * cos(θv), au * cos(θu) - av * sin(θv))

        push!(
            ellipses,
            TidalEllipse(
                eastward.constituents[k],
                wp + wm,
                wp - wm,
                mod(rad2deg((θm - θp) / 2), 180),
                mod(rad2deg(-(θp + θm) / 2), 360),
            ),
        )
    end

    return ellipses
end

"""
    quantiles(values, probabilities)

The empirical quantiles of `values`, by linear interpolation between order statistics — enough for
the QQ plots of METreport 11/2017 Figs. 21 and 30 and small enough not to want a package for.
Non-finite values are dropped.
"""
function quantiles(values, probabilities)
    sorted = sort!(Float64[value for value in values if isfinite(value)])
    isempty(sorted) && return fill(NaN, length(probabilities))

    n = length(sorted)
    return map(probabilities) do p
        position = clamp(p, 0, 1) * (n - 1) + 1
        lower = floor(Int, position)
        upper = min(lower + 1, n)
        weight = position - lower
        return (1 - weight) * sorted[lower] + weight * sorted[upper]
    end
end

"""
    direction_statistics(eastward, northward)

Speed and compass direction (degrees clockwise from north, the direction the flow is going *to*) for
a pair of velocity components — what the current roses and the directional histograms of METreport
11/2017 Figs. 19 and 20 are built from.
"""
function direction_statistics(eastward, northward)
    speed = Float64[]
    direction = Float64[]
    for (u, v) in zip(eastward, northward)
        if isfinite(u) && isfinite(v)
            push!(speed, hypot(u, v))
            push!(direction, mod(rad2deg(atan(u, v)), 360))
        else
            push!(speed, NaN)
            push!(direction, NaN)
        end
    end

    return speed, direction
end
