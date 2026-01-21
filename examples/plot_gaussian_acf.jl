# Note: Requires Plots.jl and LaTeXStrings.jl

using SommerfeldCZT
using Plots
using LaTeXStrings

# 1. Physical Parameters
f(k) = exp(-k^2 / 2.0)  # Gaussian Power Spectral Density
kmax = 8.0
x = range(0, 5, length=100) # Time-lags (τ)

# 2. Computation
# real.() because the ACF of a symmetric real spectrum is real
results = real.(sommerfeld_integral(f, x, kmax, 10))

# 3. Analytical Reference: I(x) = sqrt(π/2) * exp(-x^2/2)
analytical(x) = sqrt(π / 2) * exp(-x^2 / 2)

# 4. Plotting
plot(x, results,
    label="SommerfeldCZT",
    lw=2,
    color=:blue,
    bg=:white)

scatter!(x[1:5:end], analytical.(x[1:5:end]),
    label="Analytical Reference",
    markershape=:circle,
    color=:red)

xlabel!(L"Dimensionless Lag $x$")
ylabel!(L"ACF Amplitude $I(x)$")
title!("Thermal Ion Line: Spectrum to Autocorrelation")

savefig("gaussian_acf_validation.png")