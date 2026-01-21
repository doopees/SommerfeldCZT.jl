module SommerfeldCZT

using FFTW

export sommerfeld_integral

# Internal helper to handle the convolution steps
function _czt_convolve!(x_total, y_total, Xₙ, Yₙ, L_val, offset)
    # Use views for the part of the giant buffers we need
    xₙ = view(x_total, 1:L_val)
    yₙ = view(y_total, 1:L_val)

    # Reset views to zero
    fill!(xₙ, 0)
    fill!(yₙ, 0)

    # Pre-plan the FFTs
    p_fft = plan_fft(xₙ)
    p_ifft = plan_ifft(xₙ)

    # Fill buffers
    xₙ[1:length(Xₙ)] .= Xₙ
    yₙ[1:length(Yₙ)] .= Yₙ[end:-1:1]
    if length(Yₙ) > 1
        yₙ[length(Yₙ)+1:end] .= Yₙ[2:end]
    end

    # FFT-based high-speed convolution
    conv_full = p_ifft * ((p_fft * xₙ) .* (p_fft * yₙ))

    # Return the relevant N range samples based on the caller's required offset
    return conv_full[offset:offset+length(Xₙ)-1]
end

"""
    sommerfeld_integral(f, x, kmax, iters=1, α=-im)

Evaluates the Sommerfeld-type integral ∫ f(k) exp(α k x) dk using the 
Adaptive Chirp-Z Transform and Simpson's rule refinement as described by Li (1991).
"""
function sommerfeld_integral(f::Function, x::AbstractVector{T}, kₘₐₓ::Real, iters::Int=1, α::Complex=Complex{T}(-1.0im)) where T<:Real
    N = length(x)
    x₀ = x[1]
    Δx = x[2] - x₀
    Δk = kₘₐₓ / N

    # Precompute indices and Chirp-Z parameter W
    m = 0:N-1
    n = 0:N-1
    W = exp(α * Δk * Δx)

    # --- 1. Pre-allocation ---
    max_X_len = (iters == 1) ? N : (2^(iters - 2) * N)
    max_L = 2 * max_X_len - 1

    x_total = zeros(Complex{T}, max_L)
    y_total = zeros(Complex{T}, max_L)

    # --- 2. Initial Stage (Trapezoidal Rule) ---
    k = range(0, kₘₐₓ, length=N + 1)
    fₖ = f.(k)
    gₙ = [fₖ[1] / 2; fₖ[2:end-1]]

    Xₙ = gₙ .* exp.(α .* n .* Δk .* x₀) .* W .^ (n .^ 2 ./ 2)
    Yₙ = W .^ (-n .^ 2 ./ 2)

    # Initial convolution (offset is N)
    conv_XₙYₙ = _czt_convolve!(x_total, y_total, Xₙ, Yₙ, 2N - 1, N)

    R = Δk * (fₖ[end] / 2) .* exp.(α .* kₘₐₓ .* x)
    S = Δk .* W .^ (m .^ 2 ./ 2) .* conv_XₙYₙ .+ R

    if iters == 1
        return S
    end

    # --- 3. Adaptive Iterations (Simpson's Rule) ---
    for i in 2 .^ (1:iters-1)
        nₘₐₓ = (i ÷ 2) * N - 1
        n = 0:nₘₐₓ
        Δkᵢ = Δk / i

        gₙ = f.((2 .* n .+ 1) .* Δkᵢ)

        Xₙ = gₙ .* exp.(α .* (2 .* n .+ 1) .* Δkᵢ .* x₀) .* W .^ (n .^ 2 ./ i)
        Yₙ = W .^ (-n .^ 2 ./ i)

        # Adaptive convolution (offset is n_max + 1)
        conv_XₙYₙ = _czt_convolve!(x_total, y_total, Xₙ, Yₙ, 2 * length(Xₙ) - 1, nₘₐₓ + 1)

        # Refine S using the midpoint samples
        Sᵢ = S ./ 2 .+ Δkᵢ .* W .^ ((m .^ 2 .+ m) ./ i) .* conv_XₙYₙ[1:N]
        S = (4 .* Sᵢ .- S) ./ 3  # Simpson's Rule
    end

    return S
end

end # module