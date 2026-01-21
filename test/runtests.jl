using Test
using SommerfeldCZT
using SpecialFunctions

@testset "SommerfeldCZT.jl Comprehensive Tests" begin

    # 1. Analytical Validation: Exponential Kernel
    # Verification against a basic integral with a known analytical solution:
    # ∫_0^1 exp(-ik) exp(-im*k*x) dk = (exp(-im(1+x)) - 1) / -im(1+x)
    @testset "Exponential Kernel" begin
        f(k) = exp(-1.0im * k)
        kmax = 1.0
        x_vec = collect(range(0.1, 1.0, length=10))

        res = sommerfeld_integral(f, x_vec, kmax, 10)

        # Analytical solution
        analytical(x) = (exp(-1.0im * (1 + x)) - 1) / (-1.0im * (1 + x))

        @test res ≈ analytical.(x_vec) atol = 1e-5
    end

    # 2. Physical Significance: ISR Ion Line (Thermal Gaussian)
    # The simplest ISR model assumes a Gaussian velocity distribution.
    # Spectral Density: S(k) = exp(-k²/2)
    # ACF (Time Domain): R(τ) = exp(-τ²/2)
    # This is the "Thermal limit" of an ISR signal.
    @testset "Thermal Ion Line" begin
        f(k) = exp(-k^2 / 2.0)

        # Analytical 1D transform (one-sided integration 0 to ∞):
        # ∫₀∞ exp(-k²/2) cos(kx) dk = √(π/2) * exp(-x²/2)
        analytical_thermal(x) = sqrt(pi / 2.0) * exp(-x^2 / 2.0)

        x_vec = collect(range(0.1, 2.0, length=5))
        kmax = 8.0 # Gaussian decays so fast that kmax=8 is effectively ∞
        res = sommerfeld_integral(f, x_vec, kmax, 10)

        @test real.(res) ≈ analytical_thermal.(x_vec) atol = 1e-4
    end

    # 3. Numerical precision: Comparison with Adaptive Gauss-Kronrod
    # Stress-test the CZT logic against a high-precision Adaptive 
    # Gauss-Kronrod integrator (QuadGK) for a damped oscillatory integrand.
    using QuadGK
    @testset "Comparison with QuadGK" begin
        # A typical ISR-like integrand: damped oscillator
        f(k) = exp(-0.5 * k) * cos(k)
        kmax = 10.0
        r_vec = collect(range(0.1, 2.0, length=5))

        # Calculate using SommerfeldCZT
        res_czt = sommerfeld_integral(f, r_vec, kmax, 15)

        # Calculate using high-precision QuadGK for each r
        # We integrate: f(k) * exp(-im * k * r) dk from 0 to kmax
        function reference_quad(r)
            val, err = quadgk(k -> f(k) * exp(-1.0im * k * r), 0, kmax)
            return val
        end
        res_quad = reference_quad.(r_vec)

        @test res_czt ≈ res_quad atol = 1e-5
    end

    # 4. Mathematical Consistency: Linearity Test
    # Verify that the integral maintains linearity
    @testset "Linearity" begin
        f1(k) = exp(-k)
        f2(k) = exp(-2k)
        r = [1.0, 2.0]
        kmax = 5.0

        res1 = sommerfeld_integral(f1, r, kmax, 2)
        res2 = sommerfeld_integral(f2, r, kmax, 2)
        res_both = sommerfeld_integral(k -> f1(k) + f2(k), r, kmax, 2)

        @test res_both ≈ (res1 + res2) atol = 1e-12
    end

    # 5. Edge Cases
    # Verify behavior for edge cases like zero integrand
    @testset "Edge Cases" begin
        # Zero integrand should return zero
        @test all(sommerfeld_integral(k -> 0.0, [1.0, 5.0], 10.0, 1) .== 0.0)
    end
end