# Motor-unit accessors: discharge pulse trains, firing statistics, and spatial
# action-potential maps derived from an hdemg_decomposition.

#' Motor-unit discharge pulse trains
#'
#' Returns the discharge instants of every decomposed motor unit in tidy form.
#'
#' @param decomposition An `hdemg_decomposition` from [hdEMGDecompose()].
#' @return A data frame with one row per discharge: `unit` (integer),
#'   `discharge` (sample index), and `time_sec` (discharge time in seconds).
#' @seealso [hdEMGDecompose()], [muFiringStats()]
#' @export
#' @examples
#' sim <- make_hdemg_sim(n_units = 2, duration_sec = 3, grid = c(4, 4),
#'                       noise_sd = 0.005, seed = 1)
#' dec <- hdEMGDecompose(sim, n_units = 4, extension_factor = 10, seed = 1)
#' head(muPulseTrains(dec))
muPulseTrains <- function(decomposition) {
  stopifnot(inherits(decomposition, "hdemg_decomposition"))
  sr <- decomposition$sampling_rate
  pt <- decomposition$pulse_trains
  if (!length(pt)) {
    return(data.frame(unit = integer(0), discharge = integer(0),
                      time_sec = numeric(0)))
  }
  do.call(rbind, lapply(seq_along(pt), function(u) {
    data.frame(unit = u, discharge = pt[[u]],
               time_sec = (pt[[u]] - 1) / sr, stringsAsFactors = FALSE)
  }))
}

#' Motor-unit firing statistics
#'
#' Per-unit discharge-rate and recruitment statistics: the mean discharge rate,
#' the mean and coefficient of variation of the inter-spike interval, the
#' recruitment (first discharge) and derecruitment (last discharge) times, and
#' the decomposition quality (pulse-to-noise ratio and silhouette).
#'
#' @param decomposition An `hdemg_decomposition` from [hdEMGDecompose()].
#' @return A data frame with one row per motor unit and columns `unit`,
#'   `n_discharges`, `discharge_rate` (Hz), `mean_isi_sec`, `cov_isi`,
#'   `recruitment_sec`, `derecruitment_sec`, `pnr_db`, `sil`.
#' @references De Luca, C.J. et al. (1982). Behaviour of human motor units in
#'   different muscles during linearly varying contractions. J Physiol 329:113-128.
#' @seealso [hdEMGDecompose()], [muPulseTrains()]
#' @export
#' @examples
#' sim <- make_hdemg_sim(n_units = 2, duration_sec = 3, grid = c(4, 4),
#'                       noise_sd = 0.005, seed = 1)
#' dec <- hdEMGDecompose(sim, n_units = 4, extension_factor = 10, seed = 1)
#' muFiringStats(dec)
muFiringStats <- function(decomposition) {
  stopifnot(inherits(decomposition, "hdemg_decomposition"))
  sr <- decomposition$sampling_rate
  dur <- decomposition$n_time / sr
  pt <- decomposition$pulse_trains
  q <- decomposition$quality
  if (!length(pt)) {
    return(data.frame(
      unit = integer(0), n_discharges = integer(0), discharge_rate = numeric(0),
      mean_isi_sec = numeric(0), cov_isi = numeric(0),
      recruitment_sec = numeric(0), derecruitment_sec = numeric(0),
      pnr_db = numeric(0), sil = numeric(0)))
  }
  rows <- lapply(seq_along(pt), function(u) {
    sp <- sort(pt[[u]])
    isi <- diff(sp) / sr
    data.frame(
      unit = u, n_discharges = length(sp),
      discharge_rate = length(sp) / dur,
      mean_isi_sec = if (length(isi)) mean(isi) else NA_real_,
      cov_isi = if (length(isi) > 1L && mean(isi) > 0) {
        stats::sd(isi) / mean(isi)
      } else {
        NA_real_
      },
      recruitment_sec = (sp[1] - 1) / sr,
      derecruitment_sec = (sp[length(sp)] - 1) / sr,
      pnr_db = q$pnr_db[u], sil = q$sil[u], stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Motor-unit action-potential maps
#'
#' Spike-triggered average of the raw HD-sEMG for each motor unit, yielding the
#' multi-channel action-potential waveform and its peak-to-peak amplitude across
#' the electrode grid (the spatial "map" used to localise the unit).
#'
#' @param decomposition An `hdemg_decomposition` from [hdEMGDecompose()].
#' @param window_ms Half-window in milliseconds averaged around each discharge
#'   (default 15).
#' @return A length-`n_units` list; each element has `waveform`
#'   (`(2*w+1) x n_channels` spike-triggered average), `amplitude` (per-channel
#'   peak-to-peak), `map` (the amplitude reshaped to the `rows x cols` grid, or
#'   `NULL` if the grid is unknown), and `peak_channel` (the strongest channel).
#' @references Farina, D. et al. (2016). Reflections on the evolution of
#'   surface EMG decomposition. J Electromyogr Kinesiol 26:39-46.
#' @seealso [hdEMGDecompose()]
#' @export
#' @examples
#' sim <- make_hdemg_sim(n_units = 2, duration_sec = 3, grid = c(4, 4),
#'                       noise_sd = 0.005, seed = 1)
#' dec <- hdEMGDecompose(sim, n_units = 4, extension_factor = 10, seed = 1)
#' maps <- muActionPotentialMaps(dec)
#' maps[[1]]$peak_channel
muActionPotentialMaps <- function(decomposition, window_ms = 15) {
  stopifnot(inherits(decomposition, "hdemg_decomposition"))
  stopifnot(is.numeric(window_ms), length(window_ms) == 1L, window_ms > 0)
  sig <- decomposition$signal
  sr <- decomposition$sampling_rate
  grid <- decomposition$grid
  N <- nrow(sig); M <- ncol(sig)
  w <- max(1L, as.integer(round(window_ms / 1000 * sr)))
  offs <- (-w):w

  lapply(decomposition$pulse_trains, function(sp) {
    sp <- sp[sp > w & sp <= (N - w)]
    if (!length(sp)) {
      wave <- matrix(NA_real_, length(offs), M)
    } else {
      acc <- matrix(0, length(offs), M)
      for (s in sp) acc <- acc + sig[s + offs, , drop = FALSE]
      wave <- acc / length(sp)
    }
    amp <- apply(wave, 2, function(v) max(v) - min(v))
    map <- if (!is.null(grid) && length(grid) == 2L &&
               prod(grid) == M) {
      matrix(amp, nrow = grid[1], ncol = grid[2])
    } else {
      NULL
    }
    list(waveform = wave, amplitude = amp, map = map,
         peak_channel = if (all(is.na(amp))) NA_integer_ else which.max(amp))
  })
}
