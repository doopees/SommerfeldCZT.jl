# Note: Requires QuadGK.jl and BenchmarkTools.jl

using SommerfeldCZT
using QuadGK
using BenchmarkTools

# A harder, oscillatory integrand
f(k) = exp(-0.1 * k) * cos(10.0 * k)
kmax = 50.0  # Larger range to force more work
x_vec = collect(range(0.1, 10.0, length=100_000)) # High resolution (100k points)

println("Benchmarking CZT (100k points)...")
b_czt = @benchmark sommerfeld_integral($f, $x_vec, $kmax, 1)

println("Benchmarking QuadGK (100k points)...")
# Note: we use a lower tolerance for QuadGK to give it a fighting chance
b_quad = @benchmark [quadgk(k -> $f(k) * exp(-1.0im * k * x), 0, $kmax, rtol=1e-3)[1] for x in $x_vec]

speedup = median(b_quad.times) / median(b_czt.times)
println("Real Speedup factor: ", speedup)