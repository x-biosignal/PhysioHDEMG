# PhysioHDEMG 0.2.1

## Documentation

* `?PhysioHDEMG` now answers: a package help page gives one paragraph on what the
  package is for, the main entry points grouped by task, and where to go next.
* A vignette carries one task end to end on synthetic or bundled data, offline,
  and is built and run by `R CMD check`.
* Runnable `@examples` added or corrected across 4 help pages. Each runs
  offline in seconds, writes nothing outside `tempdir()`, and is executed by
  `R CMD check`; anything needing a device, a download or an optional backend is
  fenced with the reason stated.

# PhysioHDEMG 0.2.0

- `readHDEMG()`: read an HD-sEMG electrode-grid recording from an HDF5 file (a
  `signal` dataset and a `sampling_rate`, an optional `grid`, and an optional
  `reference_mus` group of per-unit discharge-sample vectors) into a
  `PhysioExperiment`, orienting the signal to time x channels and carrying any
  bundled reference decomposition in metadata `reference_mupulses`. The on-disk
  ingest counterpart of `hdemgToPhysioExperiment()`; requires `rhdf5`.
- `matchMotorUnits()`: match a decomposition to a reference motor-unit set with a
  lag-aware rate-of-agreement -- for each reference unit, the best-matching
  decomposed unit (allowing a constant latency up to `max_lag`) and the fraction
  of discharges the two share within `tol`. Scores cross-tool *recovery* ("which
  reference units did we recover, and which did we miss") rather than
  bit-identical agreement.


# PhysioHDEMG 0.1.1

## Validation & features

- Validated `hdEMGDecompose()` against **openhdemg** on a real 64-channel HD-sEMG
  recording (VAL-12/VAL-13): it recovers **4 of 5 reference motor units at
  rate-of-agreement 0.79-0.99** (lag-aligned). (An earlier report of ~0.02
  agreement was a measurement artifact from omitting the constant decomposition
  lag when matching discharge times; lag alignment is standard and PhysioHDEMG's
  own `test-decompose.R` already uses it.)
- Added a `bandpass` argument to `hdEMGDecompose()` for zero-phase band-pass
  conditioning of real HD-sEMG (e.g. `bandpass = c(20, 500)`). Documented that
  real recordings need lower quality thresholds (`pnr_threshold ~ 20`,
  `sil_threshold ~ 0.85`) than the clean-data defaults.

# PhysioHDEMG 0.1.0

- Initial release: HD-sEMG I/O, simulation, convolutive blind-source-separation
  motor-unit decomposition, firing statistics, and action-potential maps.
