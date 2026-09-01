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

#' Read a high-density surface EMG recording from an HDF5 file
#'
#' Reads an HD-sEMG electrode-grid recording from an HDF5 (`.h5`) file into a
#' [PhysioExperiment][hdemgToPhysioExperiment], the ingest counterpart of
#' [hdemgToPhysioExperiment()] for on-disk grids. The file must contain a
#' `signal` dataset (the raw grid signal) and a `sampling_rate` scalar; it may
#' also contain a `grid` dataset (`c(rows, cols)`) and a `reference_mus` group of
#' per-unit integer discharge-sample vectors -- a bundled reference decomposition
#' (e.g. exported from the openhdemg sample file) -- which are stored in the
#' returned object's metadata as `reference_mupulses` for scoring a decomposition
#' with [matchMotorUnits()].
#'
#' The `signal` dataset may be stored `n_time x n_channels` or
#' `n_channels x n_time`; the longer axis is taken as time.
#'
#' @param path Path to the HDF5 file.
#' @param grid Optional integer `c(rows, cols)` electrode-grid layout override;
#'   taken from the file's `grid` dataset when present, otherwise
#'   `c(1, n_channels)`.
#' @param signal_name,rate_name Dataset names for the signal and the sampling
#'   rate (defaults `"signal"` and `"sampling_rate"`).
#' @param unit Signal unit stored in the channel metadata (default `"uV"`).
#' @return A `PhysioExperiment` holding the HD-sEMG grid, with any bundled
#'   reference decomposition in metadata `reference_mupulses`.
#' @seealso [hdEMGDecompose()], [matchMotorUnits()], [hdemgToPhysioExperiment()]
#' @export
readHDEMG <- function(path, grid = NULL, signal_name = "signal",
                      rate_name = "sampling_rate", unit = "uV") {
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("readHDEMG() requires the rhdf5 package.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(sprintf("HD-EMG file not found: %s", path), call. = FALSE)
  }
  contents <- rhdf5::h5ls(path)
  sig <- as.matrix(rhdf5::h5read(path, signal_name))
  storage.mode(sig) <- "double"
  if (ncol(sig) > nrow(sig)) sig <- t(sig)          # orient to n_time x n_channels
  sr <- as.numeric(rhdf5::h5read(path, rate_name))[1]
  n_ch <- ncol(sig)
  if (is.null(grid)) {
    grid <- if ("grid" %in% contents$name) as.integer(rhdf5::h5read(path, "grid"))
            else c(1L, n_ch)
  }
  grid <- as.integer(grid)
  if (length(grid) != 2L || prod(grid) != n_ch) grid <- c(1L, n_ch)
  refmu <- NULL
  nm <- contents$name[contents$group == "/reference_mus"]
  if (length(nm)) {
    nm <- nm[order(as.integer(gsub("\\D", "", nm)))]  # numeric order (mu_1, mu_2, ...)
    refmu <- lapply(nm, function(n)
      as.integer(rhdf5::h5read(path, paste0("reference_mus/", n))))
  }
  pe <- hdemgToPhysioExperiment(sig, sampling_rate = sr, grid = grid, unit = unit)
  if (!is.null(refmu)) {
    md <- S4Vectors::metadata(pe)
    md$reference_mupulses <- refmu
    S4Vectors::metadata(pe) <- md
  }
  pe
}
