# Ingest of HD-sEMG electrode-grid recordings into PhysioExperiment objects.

#' Wrap an HD-sEMG electrode grid in a PhysioExperiment
#'
#' Builds a `PhysioExperiment` (from PhysioCore) holding a single HD-sEMG assay
#' (`n_time x n_channels`) at the given sampling rate, recording the electrode
#' grid layout (`rows`, `cols`) in the channel metadata and object metadata so
#' that [hdEMGDecompose()] can recover the grid for action-potential maps.
#'
#' @param signal An `n_time x n_channels` matrix of HD-sEMG (channels ordered
#'   column-major over the grid: row varies fastest).
#' @param sampling_rate Sampling rate in Hz.
#' @param grid Integer `c(rows, cols)` electrode-grid layout; `rows * cols` must
#'   equal the number of channels.
#' @param unit Signal unit stored in the channel metadata (default `"uV"`).
#' @param metadata Optional named list stored in the object metadata.
#' @return A `PhysioExperiment` object.
#' @seealso [hdEMGDecompose()], [make_hdemg_sim()]
#' @export
hdemgToPhysioExperiment <- function(signal, sampling_rate, grid,
                                    unit = "uV", metadata = list()) {
  signal <- as.matrix(signal)
  if (!is.numeric(signal) || nrow(signal) < 3L || ncol(signal) < 2L) {
    stop("`signal` must be a numeric n_time x n_channels matrix (>= 2 ",
         "channels).", call. = FALSE)
  }
  if (!is.numeric(sampling_rate) || length(sampling_rate) != 1L ||
      !is.finite(sampling_rate) || sampling_rate <= 0) {
    stop("`sampling_rate` must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(grid) || length(grid) != 2L || prod(grid) != ncol(signal)) {
    stop("`grid` must be c(rows, cols) with rows * cols == ncol(signal).",
         call. = FALSE)
  }
  if (!requireNamespace("PhysioCore", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("hdemgToPhysioExperiment() requires the PhysioCore and S4Vectors ",
         "packages.", call. = FALSE)
  }
  grid <- as.integer(grid)
  ecoord <- expand.grid(row = seq_len(grid[1]), col = seq_len(grid[2]))
  labels <- sprintf("E%d_%d", ecoord$row, ecoord$col)
  colnames(signal) <- labels

  PhysioCore::PhysioExperiment(
    assays = S4Vectors::SimpleList(hdemg = signal),
    colData = S4Vectors::DataFrame(label = labels, type = "hdEMG",
                                   grid_row = ecoord$row,
                                   grid_col = ecoord$col, unit = unit),
    samplingRate = sampling_rate,
    metadata = c(metadata, list(hdemg_grid = grid))
  )
}
