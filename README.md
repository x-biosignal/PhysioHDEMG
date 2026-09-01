# PhysioHDEMG

[![r-universe](https://x-biosignal.r-universe.dev/badges/PhysioHDEMG)](https://x-biosignal.r-universe.dev/PhysioHDEMG)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

`PhysioHDEMG` decomposes high-density surface EMG electrode-grid recordings
into motor-unit discharge patterns. It includes simulation, pulse-train and
firing-statistic accessors, spatial motor-unit action-potential maps, and an
optional bridge to the shared `PhysioExperiment` data model.

## Installation

```r
options(repos = c(
  xbiosignal = "https://x-biosignal.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
))
install.packages("PhysioHDEMG")
```

## Quick start

```r
library(PhysioHDEMG)

sim <- make_hdemg_sim(
  n_units = 2,
  duration_sec = 3,
  grid = c(4, 4),
  noise_sd = 0.005,
  seed = 1
)

decomposition <- hdEMGDecompose(
  sim,
  n_units = 4,
  extension_factor = 10,
  seed = 1
)

muFiringStats(decomposition)
```

## Main functions

| Function | Purpose |
|---|---|
| `make_hdemg_sim()` | Generate a grid recording with known motor units |
| `hdEMGDecompose()` | Estimate accepted motor-unit discharge patterns |
| `muPulseTrains()` | Return discharges in a tabular form |
| `muFiringStats()` | Summarize discharge rate, timing, and quality measures |
| `muActionPotentialMaps()` | Recover spatial and waveform summaries by unit |
| `hdemgToPhysioExperiment()` | Convert a matrix to the shared ecosystem model |

Decomposition quality depends on the recording, electrode layout, signal
quality, and thresholds. Returned units and quality measures should be
validated for the intended research protocol.

## Ecosystem role

The package owns HD-sEMG motor-unit decomposition. Conventional EMG features
remain in `PhysioEMG`; common data structures and provenance are supplied by
`PhysioCore` when that optional integration is used.

## Documentation

- [Function reference](https://x-biosignal.r-universe.dev/PhysioHDEMG)
- [Source repository](https://github.com/x-biosignal/PhysioHDEMG)
- [Issue tracker](https://github.com/x-biosignal/PhysioHDEMG/issues)

## Citation

```r
citation("PhysioHDEMG")
```

See the ecosystem [governance](https://github.com/x-biosignal/PhysioExperiment/blob/main/GOVERNANCE.md),
[support policy](https://github.com/x-biosignal/PhysioExperiment/blob/main/SUPPORT.md),
and [contribution guide](https://github.com/x-biosignal/PhysioExperiment/blob/main/CONTRIBUTING.md).

## Author and license

Author and maintainer: **Yusuke Matsui**. Licensed under the [MIT License](LICENSE).
