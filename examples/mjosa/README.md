# Mjøsa

Facts about the data `mjosa.jl` uses. Problems with the simulation itself are in
[`KNOWN_ISSUES.md`](KNOWN_ISSUES.md).

## Rivers

Gudbrandsdalslågen, three quarters of the inflow, has no 2025 water temperature after August
(Hovdefossen `2.656.0` ends 2025-08-19), so it carries Glomma's at Glomstadfoss. Glomma reads down
to −1.3 °C in winter; `add_rivers` raises every river temperature below 0 °C to 0 °C.
