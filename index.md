---
layout: post
title: SommerfeldCZT.jl
use_math: true
---

<script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
<script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
<script>
  window.MathJax = {
    tex: {
      inlineMath: [['$', '$'], ['\\(', '\\)']],
      displayMath: [['$$', '$$']]
    }
  };
</script>

# SommerfeldCZT.jl

[![Build Status](https://github.com/doopees/SommerfeldCZT.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/doopees/SommerfeldCZT.jl/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Julia 1.10+](https://img.shields.io/badge/julia-1.10+-cb3c33.svg)](https://julialang.org/)

A high-performance Julia implementation of the **Adaptive Chirp-Z Transform (CZT)** for evaluating Sommerfeld-type integrals, based on the method by Li et al. [[1]](#1). Achieved up to **90x faster** evaluation and **97% memory reduction** compared to standard iterative quadrature methods.

## Overview
Sommerfeld integrals are common in electromagnetics and radar science. They are often highly oscillatory and difficult to compute using standard quadrature. This package uses the Chirp-Z Transform to convert the integral into a convolution, which is then solved efficiently using FFTs.

## Mathematical Foundation

The package solves integrals of the general form

$$I(x) = \int_{0}^{k_{max}} f(k) e^{\alpha k x} dk.$$

Where $x$ is typically a spatial coordinate or time delay, and $\alpha$ is an arbitrary complex constant. By discretizing the integral into $N$ samples with a spacing of $Δk=k_{max}​/N$, the integral is approximated as a Riemann sum:

$$I(x_m) \approx Δk \sum_{n=0}^{N-1} f(nΔk) e^{\alpha nΔk x_m},$$

for $x_m$ an arbitrary  vector of samples of $x$. While a standard DFT requires $x_m$ to be tied to the reciprocal of the sampling rate, `SommerfeldCZT` uses the CZT to evaluate this sum for any linear range of $x_m$ with $O(N\log N)$ complexity. See [[1]](#1) for full details.

## Key Features
* **Computational Efficiency:** Optimized for large-scale data ($N>10^3$). Operations are reduced to $O(N\log N)$ [[1]](#1), improving on the $O(N^2)$ growth of point-by-point quadrature (up to **90x speedup** in testing).
* **Extreme Memory Efficiency:** Uses pre-allocated buffers and in-place `FFTW` transforms to reduce memory overhead compared to point-by-point quadrature (up to **97% reduction** in testing).
* **Adaptive Refinement:** Iteratively doubles the number of integration samples until convergence.
* **Contour Flexibility:** The complex parameter $\alpha$ allows the evaluation to occur along paths in the complex plane, essential for bypassing poles or branch cuts in integrands.
* **Lightweight:** Minimal dependencies (only requires `FFTW.jl`).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/doopees/SommerfeldCZT.jl")
```

## API Reference

The package provides a single function for evaluating Sommerfeld integrals.

### `sommerfeld_integral`

Computes $I(x)$ for an entire vector `x`. Returns a `Vector{Complex}` of the same length as the input `x`.

| Argument | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `f` | `Function` | — | The complex integrand $f(k)$ to be evaluated. |
| `x` | `AbstractVector` | — | The target positions in $x$-space (e.g., time or distance) for evaluation. |
| `kₘₐₓ` | `Real` | — | The upper limit of integration ($0$ to $k_{max}$). |
| `iters` | `Int` | `1` | The number of adaptive refinement steps (each step doubles the sampling density). |
| `α` | `Complex` | `-1.0im` | The complex factor in the exponent $e^{\alpha k x}$. |

## Usage

You can use the package for complex cases involving contour deformation by using the `α` parameter.

```julia
using SommerfeldCZT

# 1. Define the integrand (e.g., a thermal ion line Gaussian)
f(k) = exp(-k^2 / 2.0)

# 2. Define the target positions (e.g., radar time-delays)
r = collect(range(0.1, 2.0, length=100))

# 3. Simple evaluation
# Uses kmax=8.0 and default iters=1, α=-1.0im
results = sommerfeld_integral(f, r, 8.0)

# 4. High-precision evaluation with adaptive refinement
# Doubling sampling density 10 times
refined_results = sommerfeld_integral(f, r, 8.0, 10)

# 5. Advanced evaluation with custom propagation constant α
# Useful for shifting the integration path to avoid poles
alpha_custom = -1.0im + 0.01 
custom_results = sommerfeld_integral(f, r, 8.0, 5, alpha_custom)
```

## Numerical Validation

The `sommerfeld_integral` function evaluates the integral $I(x)$ across the entire range of observation points `x` simultaneously. To verify accuracy, we compare the CZT output against the adaptive Gauss-Kronrod method ([QuadGK](https://juliamath.github.io/QuadGK.jl/stable/)) for a damped oscillatory integrand

$$f(k) = e^{-0.5k} \cos(k).$$

The following table shows 5 representative samples from a vectorized computation using 15 refinement iterations:

| Observation Point ($x$) | CZT Result (Real Part) | QuadGK (Real Part) | Absolute Error |
|:--- |:--- |:--- |:--- |
| 0.100 | 0.40717 | 0.40717 | $1.11 \times 10^{-7}$ |
| 0.575 | 0.67141 | 0.67141 | $1.16 \times 10^{-7}$ |
| 1.050 | 1.04229 | 1.04229 | $1.31 \times 10^{-7}$ |
| 1.525 | 0.50872 | 0.50872 | $1.82 \times 10^{-7}$ |
| 2.000 | 0.22558 | 0.22558 | $3.03 \times 10^{-7}$ |

The algorithm is also validated against **Analytical Exponential Kernels** and **Thermal Ion Line** profiles to ensure physical consistency with standard Incoherent Scatter Radar (ISR) spectral models. See [examples/validate_quadgk.jl](https://github.com/doopees/SommerfeldCZT.jl/tree/main/examples/validate_quadgk.jl).

## Performance & Scaling

For high-resolution signal processing (e.g., $N = 10^5$ observation points), this package offers significant performance gains:

| Method | Points ($N$) | Median Execution Time | Memory Usage | Speedup |
| :--- | :---: | ---: | ---: | ---: |
| **QuadGK** | 100,000 | ~11.58 s | 1,019.3 MiB | 1x (Baseline) |
| **SommerfeldCZT.jl** | 100,000 | **129.45 ms** | **32.8 MiB** | **~89.4x** |

For smaller $N$, the overhead of FFTs may make CZT less efficient than direct quadrature, but as $N$ increases, the advantages become pronounced. See [examples/benchmark_scaling.jl](https://github.com/doopees/SommerfeldCZT.jl/tree/main/examples/benchmark_scaling.jl).

## Example: Thermal Ion Line
In Incoherent Scatter Radar (ISR) theory, the Sommerfeld integral transforms a **Power Spectral Density (PSD)** $f(k)$ into an **Autocorrelation Function (ACF)**.

The simplest model for an ion line assumes a Gaussian velocity distribution. For such a Gaussian PSD, $f(k) = e^{-k^2/2}$, the analytical Sommerfeld integral (from 0 to $\infty$) is

$$I(x) = \sqrt{\frac{\pi}{2}} e^{-x^2/2}.$$

The figure below demonstrates that the CZT implementation perfectly recovers the analytical ACF profile.

<p align="center">
  <img src="examples/gaussian_acf_validation.png" alt="Gaussian ACF Validation">
  <br>
  <i>Fig. 1. Validation of the CZT solver against an analytical Gaussian kernel.</i>
</p>

See [examples/plot_gaussian_acf.jl](https://github.com/doopees/SommerfeldCZT.jl/tree/main/examples/plot_gaussian_acf.jl).

## References
<a id="1">[1]</a> Y. L. Li, C. H. Liu and S. J. Franke, "Adaptive evaluation of the Sommerfeld-type integral using the chirp z-transform," in *IEEE Transactions on Antennas and Propagation*, vol. 39, no. 12, pp. 1788-1791, Dec. 1991, [doi: 10.1109/8.121603](https://doi.org/10.1109/8.121603).

## License
This project is licensed under the MIT License.
