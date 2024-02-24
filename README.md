# AstroForceModels.jl

[![CI](https://github.com/jmurphy6895/AstroForceModels.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/jmurphy6895/AstroForceModels.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/jmurphy6895/AstroForceModels.jl/branch/main/graph/badge.svg?token=47G4OLV6PD)](https://codecov.io/gh/jmurphy6895/AstroForceModels.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)][docs-stable-url]
[![](https://img.shields.io/badge/docs-dev-blue.svg)][docs-dev-url]
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)



This package contains the dominant astrodynamics forces affecting the orbtial trajectory of a satellite. Currently this package implents:
- [] Zonal Harmonics
- [] Solar Radiation Pressure
- [] Drag
- [] Third Body Gravity
- [] Relativistic
- [] Albedo
- [] Solid Tides

## Installation

```julia
julia> using Pkg
julia> Pkg.add("AstroForceModels")
```

## Documentation

For more information, see the [documentation][docs-stable-url].

[comment]: <>  UPDATE WITH OUR DOCS

[docs-dev-url]: https://juliaspace.github.io/SatelliteToolboxGravityModels.jl/dev
[docs-stable-url]: https://juliaspace.github.io/SatelliteToolboxGravityModels.jl/stable
