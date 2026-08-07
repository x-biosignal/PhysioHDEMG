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
