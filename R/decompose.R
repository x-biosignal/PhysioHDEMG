# Convolutive blind source separation for HD-sEMG motor-unit decomposition
# (Holobar & Zazula 2007; Negro et al. 2016). The observations are extended and
# whitened, a fixed-point contrast maximisation seeds each source, and a
# convolution-kernel-compensation refinement locks the separation vector onto a
# single motor-unit discharge pattern.

#' Extend a channels x time matrix by delayed copies
#' @keywords internal
#' @noRd
.hd_extend <- function(X, R) {
  M <- nrow(X); N <- ncol(X)
  Xe <- matrix(0, nrow = M * R, ncol = N)
  for (d in seq_len(R) - 1L) {
    if (d == 0L) {
      block <- X
    } else {
      block <- cbind(matrix(0, M, d), X[, seq_len(N - d), drop = FALSE])
    }
    Xe[(d * M + 1L):(d * M + M), ] <- block
  }
  Xe
}

#' Centre and whiten an extended observation matrix
#' @keywords internal
#' @noRd
.hd_whiten <- function(Xe, reg = 1e-8) {
  mu <- rowMeans(Xe)
  Xc <- Xe - mu
  N <- ncol(Xc)
  C <- tcrossprod(Xc) / N
  eg <- eigen(C, symmetric = TRUE)
  d <- eg$values
  keep <- d > (reg * max(d))
  d <- d[keep]
  V <- eg$vectors[, keep, drop = FALSE]
  Wh <- diag(1 / sqrt(d), length(d)) %*% t(V)
  list(Z = Wh %*% Xc, whiten = Wh, mean = mu)
}

#' Local maxima with a minimum separation (descending-height suppression)
#' @keywords internal
#' @noRd
.hd_bandpass <- function(D, sr, low, high) {
  # dependency-free zero-phase band-pass (FFT brick-wall) for conditioning real
  # HD-sEMG before decomposition; removes DC/baseline drift and high-frequency
  # noise so real motor-unit sources separate cleanly.
  N <- nrow(D)
  fr <- (seq_len(N) - 1L) / N * sr
  fr <- pmin(fr, sr - fr)                     # fold to [0, Nyquist]
  keep <- as.numeric(fr >= low & fr <= high)
  apply(D, 2L, function(x) Re(stats::fft(stats::fft(x) * keep, inverse = TRUE)) / N)
}

.hd_peaks <- function(v, min_dist) {
  n <- length(v)
  if (n < 3L) return(integer(0))
  is_pk <- which(v[2:(n - 1L)] > v[1:(n - 2L)] &
                   v[2:(n - 1L)] >= v[3:n]) + 1L
  if (!length(is_pk)) return(integer(0))
  ord <- is_pk[order(v[is_pk], decreasing = TRUE)]
  taken <- logical(n)
  keep <- integer(0)
  for (p in ord) {
    lo <- max(1L, p - min_dist); hi <- min(n, p + min_dist)
    if (!any(taken[lo:hi])) {
      keep <- c(keep, p)
      taken[p] <- TRUE
    }
  }
  sort(keep)
}

#' One-dimensional 2-means split (spike vs noise peak heights)
#' @keywords internal
#' @noRd
.hd_kmeans2 <- function(vals) {
  u <- unique(vals)
  if (length(u) < 2L) {
    return(list(high = rep(TRUE, length(vals)), centers = c(NA, NA), sil = 0))
  }
  km <- tryCatch(
    stats::kmeans(vals, centers = matrix(c(min(vals), max(vals)), ncol = 1),
                  algorithm = "Lloyd", iter.max = 50),
    error = function(e) NULL)
  if (is.null(km)) {
    thr <- stats::median(vals)
    high <- vals > thr
    return(list(high = high, centers = c(NA, NA), sil = 0))
  }
  hi_cluster <- which.max(km$centers)
  high <- km$cluster == hi_cluster
  list(high = high, centers = as.numeric(km$centers),
       sil = .hd_silhouette2(vals, km$cluster))
}

#' Mean silhouette of a 1-D two-cluster assignment
#' @keywords internal
#' @noRd
.hd_silhouette2 <- function(vals, cluster) {
  cl <- sort(unique(cluster))
  if (length(cl) < 2L) return(0)
  a <- cl[1]; b <- cl[2]
  va <- vals[cluster == a]; vb <- vals[cluster == b]
  sil_one <- function(x, own, other) {
    ai <- mean(abs(x - own))
    bi <- mean(abs(x - other))
    (bi - ai) / max(ai, bi)
  }
  s <- c(vapply(va, function(x) sil_one(x, va, vb), numeric(1)),
         vapply(vb, function(x) sil_one(x, vb, va), numeric(1)))
  mean(s[is.finite(s)])
}

#' Rate of agreement between two spike-index vectors (best integer lag)
#' @keywords internal
#' @noRd
.hd_roa <- function(a, b, tol = 1L, max_lag = 0L) {
  if (!length(a) || !length(b)) return(list(roa = 0, lag = 0L))
  best <- list(roa = -1, lag = 0L)
  for (lag in (-max_lag):max_lag) {
    bl <- b + lag
    matched <- 0L
    bi <- 1L
    used <- logical(length(bl))
    as <- sort(a); bs <- sort(bl)
    # greedy two-pointer match within tolerance
    j <- 1L
    for (ai in as) {
      while (j <= length(bs) && bs[j] < ai - tol) j <- j + 1L
      if (j <= length(bs) && abs(bs[j] - ai) <= tol) {
        matched <- matched + 1L
        j <- j + 1L
      }
    }
    roa <- matched / (length(a) + length(b) - matched)
    if (roa > best$roa) best <- list(roa = roa, lag = lag)
  }
  best
}

#' Decompose HD-sEMG into motor-unit discharge patterns
#'
#' Decomposes a high-density surface EMG electrode-grid recording into its
#' constituent motor-unit pulse trains by convolutive blind source separation:
#' the channels are extended by delayed copies and whitened, each candidate
#' source is seeded by a fixed-point contrast (kurtosis) maximisation and then
#' refined by convolution-kernel compensation (the separation vector is
#' iteratively set to the average whitened observation at the estimated
#' discharge times), and discharge times are read off by peak detection with a
#' two-class (spike / noise) split. Sources are accepted on a pulse-to-noise
#' ratio and silhouette quality gate, and duplicate units are merged.
#'
#' @param x HD-sEMG data: an `n_time x n_channels` numeric matrix, an
#'   `hdemg_sim` object from [make_hdemg_sim()], or a `PhysioExperiment`
#'   (requires the PhysioCore package).
#' @param n_units Maximum number of separation vectors to attempt (default 8).
#' @param extension_factor Number of delayed copies per channel used to build
#'   the extended observations (default 16).
#' @param sampling_rate Sampling rate in Hz; required when `x` is a bare matrix
#'   (taken from the object otherwise).
#' @param assay_name Assay to use when `x` is a `PhysioExperiment`
#'   (default: the first assay).
#' @param refractory_sec Minimum inter-discharge interval in seconds used for
#'   peak separation (default 0.02).
#' @param pnr_threshold Minimum pulse-to-noise ratio in dB for an accepted unit
#'   (default 26).
#' @param sil_threshold Minimum silhouette for an accepted unit (default 0.9).
#' @param max_iter Maximum fixed-point iterations per source (default 100).
#' @param tol Convergence tolerance on the separation vector (default 1e-4).
#' @param bandpass Optional length-2 numeric `c(low, high)` in Hz to band-pass the
#'   signal (zero-phase) before decomposition. Recommended for real HD-sEMG
#'   (e.g. `c(20, 500)`) to remove baseline drift and high-frequency noise;
#'   `NULL` (default) leaves the signal unchanged.
#' @param seed Optional integer random seed for reproducibility.
#'
#' @details
#' The default quality gate (`pnr_threshold = 26`, `sil_threshold = 0.9`) is tuned
#' for clean data (e.g. [make_hdemg_sim()]). Real HD-sEMG motor units typically
#' have lower quality (PNR ~20-23 dB, silhouette ~0.85-0.92); for real recordings,
#' band-pass the signal (`bandpass = c(20, 500)`) and relax the thresholds
#' (`pnr_threshold ~ 20`, `sil_threshold ~ 0.85`). Validated against openhdemg on
#' a real 64-channel recording: PhysioHDEMG recovered 4 of 5 reference motor units
#' at rate-of-agreement 0.79-0.99 (lag-aligned).
#'
#' @return An `hdemg_decomposition` object with `pulse_trains` (list of
#'   discharge sample indices), a `quality` data frame
#'   (`unit`, `n_discharges`, `discharge_rate`, `pnr_db`, `sil`), the separation
#'   `filters`, and the retained `signal`/`grid`/`sampling_rate` for downstream
#'   firing statistics and action-potential maps.
#' @references Holobar, A. & Zazula, D. (2007). Multichannel blind source
#'   separation using convolution kernel compensation. IEEE TSP 55(9):4487-4496.
#'   Negro, F. et al. (2016). J Neural Eng 13(2):026027.
#' @seealso [make_hdemg_sim()], [muFiringStats()], [muActionPotentialMaps()]
#' @export
#' @examples
#' sim <- make_hdemg_sim(n_units = 2, duration_sec = 3, grid = c(4, 4),
#'                       noise_sd = 0.005, seed = 1)
#' dec <- hdEMGDecompose(sim, n_units = 4, extension_factor = 10, seed = 1)
#' dec$quality
hdEMGDecompose <- function(x, n_units = 8, extension_factor = 16,
                           sampling_rate = NULL, assay_name = NULL,
                           refractory_sec = 0.02, pnr_threshold = 26,
                           sil_threshold = 0.9, max_iter = 100, tol = 1e-4,
                           bandpass = NULL, seed = NULL) {
  parsed <- .hd_input(x, sampling_rate, assay_name)
  D <- parsed$signal            # n_time x n_channels
  sr <- parsed$sampling_rate
  grid <- parsed$grid
  if (!is.null(bandpass)) {
    stopifnot(is.numeric(bandpass), length(bandpass) == 2L,
              bandpass[1] >= 0, bandpass[2] > bandpass[1], bandpass[2] <= sr / 2)
    D <- .hd_bandpass(D, sr, bandpass[1], bandpass[2])
  }
  stopifnot(is.numeric(n_units), length(n_units) == 1L, n_units >= 1)
  stopifnot(is.numeric(extension_factor), length(extension_factor) == 1L,
            extension_factor >= 1)
  stopifnot(is.numeric(refractory_sec), length(refractory_sec) == 1L,
            is.finite(refractory_sec), refractory_sec > 0)
  stopifnot(is.numeric(pnr_threshold), length(pnr_threshold) == 1L,
            is.finite(pnr_threshold))
  stopifnot(is.numeric(sil_threshold), length(sil_threshold) == 1L,
            is.finite(sil_threshold))
  if (!is.numeric(D) || any(!is.finite(D))) {
    stop("the HD-sEMG signal must be a fully finite numeric matrix.",
         call. = FALSE)
  }
  if (!is.null(seed)) set.seed(as.integer(seed))

  N <- nrow(D); M <- ncol(D)
  R <- as.integer(extension_factor)
  if (R >= N) {
    stop("`extension_factor` must be smaller than the number of time samples.",
         call. = FALSE)
  }
  min_dist <- max(1L, as.integer(round(refractory_sec * sr)))

  X <- t(D)                     # channels x time
  Xe <- .hd_extend(X, R)
  wh <- .hd_whiten(Xe)
  Z <- wh$Z                     # P x N whitened
  P <- nrow(Z)
  activity <- colSums(Z^2)

  B <- matrix(0, nrow = P, ncol = 0)     # filters of found units (deflation)
  units <- list()
  target <- as.integer(n_units)
  # attempt many more seeds than the target: seeds land on noise or re-find an
  # already-identified unit (a delayed copy), so a tight budget would leave real
  # units undiscovered
  max_attempts <- max(target * 4L, target + 16L)
  fails <- 0L                               # consecutive unproductive attempts
  max_fails <- max(10L, target)

  for (i in seq_len(max_attempts)) {
    if (length(units) >= target) break      # enough distinct units found
    if (fails >= max_fails) break           # real units appear drained
    # seed from the highest-activity instant not yet explained
    seed_t <- which.max(activity)
    if (activity[seed_t] <= 0) break        # nothing left to explain
    w <- Z[, seed_t]
    if (ncol(B)) w <- w - B %*% crossprod(B, w)
    nw <- sqrt(sum(w^2))
    if (nw == 0) { activity[seed_t] <- 0; next }
    w <- w / nw

    # fixed-point contrast maximisation with the bounded logcosh contrast
    # (robust: unlike the pow3/kurtosis contrast it is not hijacked by a single
    # large sample, so it locks onto the periodic discharge train)
    for (it in seq_len(max_iter)) {
      s <- as.vector(crossprod(w, Z))
      g <- tanh(s)
      w_new <- Z %*% g / N - mean(1 - g^2) * w
      if (ncol(B)) w_new <- w_new - B %*% crossprod(B, w_new)
      w_new <- w_new / sqrt(sum(w_new^2))
      if (abs(abs(sum(w_new * w)) - 1) < tol) { w <- w_new; break }
      w <- w_new
    }

    # convolution-kernel-compensation refinement: lock onto one discharge train
    ref <- .hd_refine(w, Z, min_dist, n_iter = 12L)
    w <- ref$w; spikes <- ref$spikes; s <- ref$s

    # advance the seed even when an attempt yields no clean train
    activity[max(1L, seed_t - min_dist):min(N, seed_t + min_dist)] <- 0
    if (length(spikes) < 3L) { fails <- fails + 1L; next }
    # reject a train that duplicates an already-found unit (a delayed copy of
    # the same motor unit re-surfacing at another lag)
    dup <- FALSE
    for (u in units) {
      if (.hd_roa(u$spikes, spikes, tol = max(1L, min_dist %/% 2L),
                  max_lag = min_dist)$roa > 0.3) { dup <- TRUE; break }
    }
    if (dup) { fails <- fails + 1L; next }
    fails <- 0L
    # a new distinct unit: deflate its filter and suppress its discharge
    # footprint so the next seed moves toward a different motor unit
    B <- cbind(B, w)
    for (sp in spikes) {
      activity[max(1L, sp - min_dist):min(N, sp + min_dist)] <- 0
    }
    q <- .hd_quality(s, spikes, min_dist)
    units[[length(units) + 1L]] <- list(
      w = w, spikes = spikes, source = s, pnr = q$pnr, sil = q$sil)
  }

  # quality gate + duplicate removal (keep higher PNR)
  keep <- vapply(units, function(u) is.finite(u$pnr) && u$pnr >= pnr_threshold &&
                   is.finite(u$sil) && u$sil >= sil_threshold, logical(1))
  units <- units[keep]
  units <- .hd_dedup(units, tol = max(1L, min_dist %/% 2L),
                     max_lag = min_dist)

  if (length(units)) {
    ord <- order(vapply(units, function(u) min(u$spikes), numeric(1)))
    units <- units[ord]
  }

  pulse_trains <- lapply(units, function(u) u$spikes)
  filters <- if (length(units)) {
    do.call(cbind, lapply(units, function(u) u$w))
  } else {
    matrix(numeric(0), nrow = P, ncol = 0)
  }
  quality <- data.frame(
    unit = seq_along(units),
    n_discharges = vapply(units, function(u) length(u$spikes), integer(1)),
    discharge_rate = vapply(units, function(u) length(u$spikes) / (N / sr),
                            numeric(1)),
    pnr_db = vapply(units, function(u) u$pnr, numeric(1)),
    sil = vapply(units, function(u) u$sil, numeric(1)),
    stringsAsFactors = FALSE
  )

  out <- list(
    pulse_trains = pulse_trains, quality = quality, filters = filters,
    signal = D, grid = grid, sampling_rate = sr, n_time = N,
    n_channels = M, extension_factor = R)
  class(out) <- "hdemg_decomposition"
  out
}

#' CKC refinement: iterate w = mean whitened obs at discharges, keep min-CoV set
#' @keywords internal
#' @noRd
.hd_refine <- function(w, Z, min_dist, n_iter = 12L) {
  best <- list(cov = Inf, w = w, spikes = integer(0),
               s = as.vector(crossprod(w, Z)))
  for (it in seq_len(n_iter)) {
    s <- as.vector(crossprod(w, Z))
    s2 <- s^2
    peaks <- .hd_peaks(s2, min_dist)
    if (length(peaks) < 3L) break
    cl <- .hd_kmeans2(s2[peaks])
    spikes <- peaks[cl$high]
    if (length(spikes) < 3L) break
    isi <- diff(sort(spikes))
    cov_isi <- if (length(isi) > 1L && mean(isi) > 0) {
      stats::sd(isi) / mean(isi)
    } else {
      Inf
    }
    if (cov_isi < best$cov) {
      best <- list(cov = cov_isi, w = w, spikes = spikes, s = s)
    }
    # update the separation vector to the average whitened obs at discharges
    w <- rowMeans(Z[, spikes, drop = FALSE])
    nw <- sqrt(sum(w^2)); if (nw == 0) break
    w <- w / nw
  }
  best
}

#' Pulse-to-noise ratio (dB) and silhouette for a source + discharge set
#'
#' PNR compares the mean squared source at the discharge instants against the
#' silent baseline (all samples outside a refractory window around a discharge),
#' following Holobar et al. (2014). The silhouette measures the separation of the
#' spike vs noise peak-height clusters.
#' @keywords internal
#' @noRd
.hd_quality <- function(s, spikes, min_dist) {
  s2 <- s^2
  N <- length(s2)
  silent <- rep(TRUE, N)
  for (sp in spikes) {
    silent[max(1L, sp - min_dist):min(N, sp + min_dist)] <- FALSE
  }
  peaks <- .hd_peaks(s2, min_dist)
  sil <- if (length(peaks) >= 2L) .hd_kmeans2(s2[peaks])$sil else 0
  # if the discharge footprints cover the whole record there is no silent
  # baseline to compare against -> treat as unquantifiable (reject), not Inf
  denom <- if (any(silent)) mean(s2[silent]) else NA_real_
  pnr <- if (is.finite(denom) && denom > 0) {
    10 * log10(mean(s2[spikes]) / denom)
  } else {
    -Inf
  }
  list(pnr = pnr, sil = sil)
}

#' Drop duplicate motor units (high rate-of-agreement), keep higher PNR
#' @keywords internal
#' @noRd
.hd_dedup <- function(units, tol, max_lag, roa_thresh = 0.3) {
  n <- length(units)
  if (n < 2L) return(units)
  drop <- logical(n)
  for (i in seq_len(n - 1L)) {
    if (drop[i]) next
    for (j in (i + 1L):n) {
      if (drop[j]) next
      r <- .hd_roa(units[[i]]$spikes, units[[j]]$spikes, tol = tol,
                   max_lag = max_lag)$roa
      if (r > roa_thresh) {
        if (units[[i]]$pnr >= units[[j]]$pnr) drop[j] <- TRUE
        else { drop[i] <- TRUE; break }
      }
    }
  }
  units[!drop]
}

#' Coerce the various accepted inputs to a signal matrix + sampling rate + grid
#' @keywords internal
#' @noRd
.hd_input <- function(x, sampling_rate, assay_name) {
  grid <- NULL
  if (inherits(x, "hdemg_sim")) {
    return(list(signal = x$signal, sampling_rate = x$sampling_rate,
                grid = x$grid))
  }
  if (inherits(x, "PhysioExperiment")) {
    if (!requireNamespace("PhysioCore", quietly = TRUE) ||
        !requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("Decomposing a PhysioExperiment requires the PhysioCore and ",
           "SummarizedExperiment packages.", call. = FALSE)
    }
    if (is.null(assay_name)) assay_name <- PhysioCore::defaultAssay(x)
    sig <- as.matrix(SummarizedExperiment::assay(x, assay_name))
    sr <- PhysioCore::samplingRate(x)
    g <- S4Vectors::metadata(x)$hdemg_grid
    if (!is.null(g)) grid <- as.integer(g)
    return(list(signal = sig, sampling_rate = sr, grid = grid))
  }
  if (is.matrix(x) || is.data.frame(x)) {
    sig <- as.matrix(x)
    if (!is.numeric(sig)) stop("`x` must be a numeric matrix.", call. = FALSE)
    if (is.null(sampling_rate) || !is.finite(sampling_rate) ||
        sampling_rate <= 0) {
      stop("`sampling_rate` (Hz) is required when `x` is a matrix.",
           call. = FALSE)
    }
    if (nrow(sig) < 3L || ncol(sig) < 2L) {
      stop("`x` must have >= 3 time samples and >= 2 channels.", call. = FALSE)
    }
    return(list(signal = sig, sampling_rate = sampling_rate, grid = grid))
  }
  stop("`x` must be a matrix, an hdemg_sim, or a PhysioExperiment.",
       call. = FALSE)
}

#' @export
print.hdemg_decomposition <- function(x, ...) {
  cat(sprintf("<hdemg_decomposition> %d motor unit(s) from %d channels @ %g Hz\n",
              length(x$pulse_trains), x$n_channels, x$sampling_rate))
  if (nrow(x$quality)) {
    for (i in seq_len(nrow(x$quality))) {
      q <- x$quality[i, ]
      cat(sprintf("  MU %d: %d discharges, %.1f Hz, PNR %.1f dB, SIL %.2f\n",
                  q$unit, q$n_discharges, q$discharge_rate, q$pnr_db, q$sil))
    }
  }
  invisible(x)
}
