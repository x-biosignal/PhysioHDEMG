# Ground-truth HD-sEMG simulator for validating motor-unit decomposition.

#' Biphasic/triphasic motor-unit action-potential template
#'
#' A unit-peak analytic MUAP waveform (scaled derivative of a Gaussian) of a
#' given length in samples.
#' @keywords internal
#' @noRd
.hd_muap_template <- function(len, width, order = 1L) {
  u <- (seq_len(len) - (len + 1) / 2) / width
  g <- exp(-u^2 / 2)
  w <- if (order == 1L) -u * g else (u^2 - 1) * g   # bi- or tri-phasic
  m <- max(abs(w))
  if (m > 0) w / m else w
}

#' Simulate a high-density surface EMG recording with known motor units
#'
#' Generates a synthetic HD-sEMG electrode-grid signal as the linear
#' superposition of a known set of motor units, each firing at a target
#' discharge rate with inter-spike-interval jitter and each with its own
#' spatial action-potential map over the grid. The ground-truth discharge
#' times are returned so that a decomposition can be validated (rate of
#' agreement, discharge-rate error).
#'
#' @param n_units Number of motor units to simulate (default 4).
#' @param duration_sec Recording duration in seconds (default 10).
#' @param sampling_rate Sampling rate in Hz (default 2048).
#' @param grid Integer `c(rows, cols)` electrode-grid layout (default
#'   `c(5, 5)`); `rows * cols` channels are generated.
#' @param firing_rates Optional numeric vector of per-unit mean discharge
#'   rates (Hz); defaults to an evenly spaced spread in `[8, 18]` Hz.
#' @param muap_length MUAP template length in samples (default 24).
#' @param cov_isi Coefficient of variation of the inter-spike interval
#'   (firing-time jitter, default 0.1).
#' @param noise_sd Standard deviation of additive white noise, in the same
#'   units as the unit-peak MUAPs (default 0.02; smaller = higher SNR).
#' @param seed Optional integer random seed for reproducibility.
#'
#' @return An `hdemg_sim` object (a list) with:
#'   \describe{
#'     \item{signal}{An `n_time x n_channels` matrix of simulated HD-sEMG.}
#'     \item{sampling_rate}{Sampling rate in Hz.}
#'     \item{grid}{The `c(rows, cols)` grid layout.}
#'     \item{spikes}{A length-`n_units` list of integer ground-truth discharge
#'       sample indices.}
#'     \item{discharge_rates}{The realised per-unit mean discharge rate (Hz).}
#'     \item{muap}{A length-`n_units` list of `muap_length x n_channels`
#'       ground-truth action-potential templates.}
#'     \item{n_units}{Number of motor units.}
#'   }
#' @references Farina, D. et al. (2010). Decoding the neural drive to muscles
#'   from the surface electromyogram. J Physiol 588(11):1969-1975.
#' @seealso [hdEMGDecompose()]
#' @export
#' @examples
#' sim <- make_hdemg_sim(n_units = 3, duration_sec = 2, grid = c(4, 4),
#'                       seed = 1)
#' dim(sim$signal)
#' lengths(sim$spikes)
make_hdemg_sim <- function(n_units = 4, duration_sec = 10, sampling_rate = 2048,
                           grid = c(5, 5), firing_rates = NULL,
                           muap_length = 24, cov_isi = 0.1, noise_sd = 0.02,
                           seed = NULL) {
  stopifnot(is.numeric(n_units), length(n_units) == 1L, n_units >= 1)
  stopifnot(is.numeric(duration_sec), length(duration_sec) == 1L,
            duration_sec > 0)
  stopifnot(is.numeric(sampling_rate), length(sampling_rate) == 1L,
            sampling_rate > 0)
  stopifnot(is.numeric(grid), length(grid) == 2L, all(grid >= 1))
  if (!is.null(seed)) set.seed(as.integer(seed))
  n_units <- as.integer(n_units)
  L <- as.integer(muap_length)
  n_rows <- as.integer(grid[1]); n_cols <- as.integer(grid[2])
  n_ch <- n_rows * n_cols
  n_time <- as.integer(round(duration_sec * sampling_rate))

  if (is.null(firing_rates)) {
    firing_rates <- if (n_units == 1L) 12 else seq(8, 18, length.out = n_units)
  }
  firing_rates <- rep_len(as.numeric(firing_rates), n_units)

  # electrode coordinates on the grid (row, col)
  ecoord <- expand.grid(row = seq_len(n_rows), col = seq_len(n_cols))

  # spread the motor-unit territories over the grid (a space-filling sub-grid
  # with small jitter) so their action-potential maps are distinguishable
  side <- ceiling(sqrt(n_units))
  cgrid <- expand.grid(
    r = seq(1, n_rows, length.out = min(side, n_rows)),
    c = seq(1, n_cols, length.out = min(side, n_cols)))
  cgrid <- cgrid[rep_len(seq_len(nrow(cgrid)), n_units), , drop = FALSE]
  centers_r <- cgrid$r + stats::runif(n_units, -0.4, 0.4)
  centers_c <- cgrid$c + stats::runif(n_units, -0.4, 0.4)

  signal <- matrix(0, nrow = n_time, ncol = n_ch)
  spikes <- vector("list", n_units)
  rates <- numeric(n_units)
  muaps <- vector("list", n_units)

  for (k in seq_len(n_units)) {
    # spatial map: Gaussian bump centred at this unit's grid territory
    ctr_r <- centers_r[k]
    ctr_c <- centers_c[k]
    sigma <- max(0.8, mean(c(n_rows, n_cols)) / 4)
    d2 <- (ecoord$row - ctr_r)^2 + (ecoord$col - ctr_c)^2
    spatial <- exp(-d2 / (2 * sigma^2))
    # temporal template: unit-specific width and phase order
    width <- stats::runif(1, L / 8, L / 5)
    tmpl <- .hd_muap_template(L, width, order = if (k %% 2L == 0L) 2L else 1L)
    amp <- stats::runif(1, 0.8, 1.2)
    muap <- outer(tmpl, spatial) * amp          # L x n_ch template
    muaps[[k]] <- muap

    # discharge times: gamma-like ISI with target rate and CoV jitter
    mean_isi <- sampling_rate / firing_rates[k]
    sd_isi <- cov_isi * mean_isi
    n_est <- as.integer(ceiling(firing_rates[k] * duration_sec * 1.4)) + 5L
    isi <- stats::rnorm(n_est, mean_isi, sd_isi)
    isi[isi < 0.5 * mean_isi] <- 0.5 * mean_isi         # refractory floor
    t0 <- stats::runif(1, L, 2 * mean_isi)
    times <- round(t0 + cumsum(isi))
    times <- times[times >= (L + 1L) & times <= (n_time - L)]
    spikes[[k]] <- as.integer(times)
    rates[k] <- length(times) / duration_sec

    # superimpose the MUAP at each discharge time (template centred on spike)
    half <- (L - 1L) %/% 2L
    for (ti in times) {
      idx <- (ti - half):(ti - half + L - 1L)
      keep <- idx >= 1L & idx <= n_time
      signal[idx[keep], ] <- signal[idx[keep], ] + muap[keep, ]
    }
  }

  if (noise_sd > 0) {
    signal <- signal + matrix(stats::rnorm(n_time * n_ch, 0, noise_sd),
                              n_time, n_ch)
  }

  out <- list(signal = signal, sampling_rate = sampling_rate,
              grid = c(n_rows, n_cols), spikes = spikes,
              discharge_rates = rates, muap = muaps, n_units = n_units)
  class(out) <- "hdemg_sim"
  out
}

#' @export
print.hdemg_sim <- function(x, ...) {
  cat(sprintf("<hdemg_sim> %d units, %d x %d grid, %.1f s @ %g Hz\n",
              x$n_units, x$grid[1], x$grid[2],
              nrow(x$signal) / x$sampling_rate, x$sampling_rate))
  cat(sprintf("  discharge rates: %s Hz\n",
              paste(sprintf("%.1f", x$discharge_rates), collapse = ", ")))
  invisible(x)
}
