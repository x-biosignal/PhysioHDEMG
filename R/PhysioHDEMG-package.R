#' PhysioHDEMG: high-density surface EMG motor-unit decomposition
#'
#' PhysioHDEMG decomposes high-density surface EMG (HD-sEMG) electrode-grid
#' recordings into their constituent motor-unit discharge patterns by
#' convolutive blind source separation, and summarises the recovered units by
#' discharge rate, recruitment, decomposition quality, and spatial
#' action-potential map over the grid. A ground-truth simulator makes the whole
#' pipeline runnable offline for validation.
#'
#' @section Simulate ground truth:
#' [make_hdemg_sim()] generates a synthetic grid recording with known motor
#' units (discharge times and spatial maps) for validating a decomposition.
#'
#' @section Decompose:
#' [hdEMGDecompose()] estimates accepted motor-unit pulse trains from a matrix,
#' an `hdemg_sim`, or a `PhysioExperiment`.
#'
#' @section Inspect motor units:
#' [muPulseTrains()] returns the discharges in tidy form, [muFiringStats()]
#' summarises discharge rate, recruitment, and quality, and
#' [muActionPotentialMaps()] recovers the per-unit spatial action-potential map.
#'
#' @section Validate against a reference:
#' [matchMotorUnits()] scores a decomposition against a reference set of
#' discharge patterns with a lag-aware rate of agreement.
#'
#' @section Ingest and interchange:
#' [readHDEMG()] reads an HD-sEMG grid (and any bundled reference decomposition)
#' from an HDF5 file, and [hdemgToPhysioExperiment()] wraps a matrix in the
#' shared ecosystem data model.
#'
#' @section Where to go next:
#' Conventional EMG features (envelopes, amplitude, spectral indices) live in
#' the sibling `PhysioEMG` package; the shared `PhysioExperiment` data model and
#' provenance come from `PhysioCore`. See `vignette("hdemg-decomposition",
#' package = "PhysioHDEMG")` for an end-to-end walk-through.
#'
#' @keywords internal
"_PACKAGE"
