# Note: Requires QuadGK.jl and QuadGK.jl

using SommerfeldCZT
using QuadGK
using Printf

# 1. Define a damped oscillatory integrand
# This is a classic test for Sommerfeld solvers
f(k) = exp(-0.5 * k) * cos(k)

kmax = 10.0
x_vec = collect(range(0.1, 2.0, length=5))

# 2. Compute using CZT (15 refinement iterations)
res_czt = sommerfeld_integral(f, x_vec, kmax, 15)

# 3. Compute using QuadGK (High-precision reference)
# We use a list comprehension to evaluate at each point x
res_quad = [quadgk(k -> f(k) * exp(-1.0im * k * x), 0, kmax)[1] for x in x_vec]

# 4. Print Comparison Table
println("-"^65)
@printf("%-10s | %-15s | %-15s | %-10s\n", "x", "CZT (Real)", "QuadGK (Real)", "Abs Error")
println("-"^65)

# eachindex(x_vec) is preferred over 1:length(x_vec)
for i in eachindex(x_vec, res_czt, res_quad)
    czt_val = real(res_czt[i])
    quad_val = real(res_quad[i])
    err = abs(czt_val - quad_val)

    @printf("%-10.3f | %-15.5f | %-15.5f | %-10.2e\n",
        x_vec[i], czt_val, quad_val, err)
end
println("-"^65)