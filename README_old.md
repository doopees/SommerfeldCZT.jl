# SommerfeldCZT.jl

A high-performance Julia implementation of the **Adaptive Chirp-Z Transform (CZT)** for evaluating Sommerfeld-type integrals, based on the method by **Li (1991)**.

## Overview
Sommerfeld integrals are common in electromagnetics and radar science (ISR). They are often highly oscillatory and difficult to compute using standard quadrature. This package uses the Chirp-Z Transform to convert the integral into a convolution, which is then solved efficiently using FFTs.

### Key Features
* **Adaptive Refinement:** Uses Simpson's rule refinement to double sampling points until convergence.
* **Optimized:** Uses pre-allocated buffers and FFTW plans to minimize overhead during iterations.
* **Unified API:** Simple `sommerfeld_integral` function for complex integrand evaluation.

## Installation
```julia
using Pkg
Pkg.add(url="[https://github.com/doopees/SommerfeldCZT.jl](https://github.com/doopees/SommerfeldCZT.jl)")