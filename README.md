[![Build Status](https://github.com/WGUNDERWOOD/DyadicKDE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/WGUNDERWOOD/DyadicKDE.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/WGUNDERWOOD/DyadicKDE.jl/branch/main/graph/badge.svg?token=JLN7VK2LT8)](https://codecov.io/gh/WGUNDERWOOD/DyadicKDE.jl)

# DyadicKDE.jl

Dyadic kernel density estimation in Julia.

## Introduction

This repository provides a Julia package which implements the methods for
dyadic kernel density estimation detailed in
[Cattaneo, Feng and Underwood, 2022](https://arxiv.org/abs/2201.05967).
In particular, the package provides the capability for computing:

- Point estimates of a dyadic density function
- Pointwise confidence intervals for the density
- Uniform confidence bands for the density
- A rule-of-thumb bandwidth for the estimator

The currently supported kernels are

- Epanechnikov, order 2
- Epanechnikov, order 4

## Installation

Install from the Julia General registry by starting a
Julia interactive session and running

```julia
] add DyadicKDE
```

Alternatively install from source by running

```julia
] add "https://github.com/WGUNDERWOOD/DyadicKDE.jl.git"
```

The package can then be loaded with

```julia
using DyadicKDE
```

and tested (this may take a few minutes) with

```julia
] test DyadicKDE
```

## Dependencies

DyadicKDE.jl depends on several other Julia packages which are listed in
[Project.toml](https://github.com/WGUNDERWOOD/DyadicKDE.jl/tree/main/Project.toml).

## Quick start guide

```julia
# load package
using DyadicKDE

# specify parameters
n_data = 100
kernel_name = "epanechnikov_order_2"
evals = collect(range(-2.0, stop=2.0, length=10))
sdp_solver = "cosmo"
n_resample = 1000
significance_level = 0.05
p = [0.25, 0.0, 0.75]

# make data and get bandwidth
data = make_data(n_data, p)
h_ROT = estimate_ROT_bandwidth(data, "epanechnikov_order_2")

# fit dyadic kernel density estimator
est = DyadicKernelDensityEstimator(
    kernel_name, h_ROT, significance_level,
    n_resample, sdp_solver, evals, data, Dict())

fit(est)

# display properties of estimator
display(est)

# display evaluation points
display(evals')

# display point estimates
display(est.fhat')

# display pointwise confidence intervals
display(est.pci)

# display uniform confidence band
display(est.ucb)
```

## Paper replication

The table of results presented in
[Cattaneo, Feng and Underwood, 2022](https://arxiv.org/abs/2201.05967)
can be replicated by running the Julia script
[test/paper_replication.jl](https://github.com/WGUNDERWOOD/DyadicKDE.jl/tree/main/test/paper_replication.jl).
This script may take a long time to run (several hours),
but can be accelerated for example by first starting Julia with
`julia -t 8`
to use 8 CPU threads.
