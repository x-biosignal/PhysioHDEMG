library(testthat)
library(PhysioHDEMG)

# self-contained rate of agreement (matched discharges within tol, best lag)
roa_lag <- function(a, b, tol, max_lag) {
  best <- 0
  as <- sort(a)
  for (lag in (-max_lag):max_lag) {
    bl <- sort(b + lag); j <- 1L; m <- 0L
    for (ai in as) {
      while (j <= length(bl) && bl[j] < ai - tol) j <- j + 1L
      if (j <= length(bl) && abs(bl[j] - ai) <= tol) { m <- m + 1L; j <- j + 1L }
    }
    r <- m / (length(a) + length(b) - m)
    if (r > best) best <- r
  }
  best
}
best_match <- function(sim, spikes, tol) {
  which.max(vapply(sim$spikes, function(s) roa_lag(s, spikes, tol, 40L),
                   numeric(1)))
}

# One high-SNR decomposition, reused across tests (decomposition is the slow
# part). noise_sd = 0.003 is the high-SNR case of the acceptance criteria.
sim <- make_hdemg_sim(n_units = 3, duration_sec = 10, grid = c(5, 5),
                      noise_sd = 0.003, seed = 1)
dec <- hdEMGDecompose(sim, n_units = 8, extension_factor = 16, seed = 1)
sr <- sim$sampling_rate
tol <- max(1L, round(0.0005 * sr))          # +/- 0.5 ms

test_that("decomposition recovers known motor units at high rate of agreement", {
  expect_s3_class(dec, "hdemg_decomposition")
  expect_gte(nrow(dec$quality), 2)          # recover most of the 3 units
  for (i in seq_along(dec$pulse_trains)) {
    r <- max(vapply(sim$spikes, function(s)
      roa_lag(s, dec$pulse_trains[[i]], tol, 40L), numeric(1)))
    expect_gt(r, 0.9)                        # > 90% rate of agreement
  }
  # recovered units are distinct (no two match the same ground-truth unit)
  matched <- vapply(dec$pulse_trains, function(p) best_match(sim, p, tol),
                    integer(1))
  expect_equal(length(unique(matched)), length(matched))
})

test_that("accepted motor units clear the PNR and silhouette gates", {
  expect_true(all(dec$quality$pnr_db > 26))
  expect_true(all(dec$quality$sil > 0.9))
})

test_that("discharge-rate estimates are within 5% of ground truth", {
  for (i in seq_along(dec$pulse_trains)) {
    j <- best_match(sim, dec$pulse_trains[[i]], tol)
    err <- abs(dec$quality$discharge_rate[i] - sim$discharge_rates[j]) /
      sim$discharge_rates[j]
    expect_lt(err, 0.05)
  }
})

test_that("muPulseTrains and muFiringStats are consistent", {
  pt <- muPulseTrains(dec)
  expect_s3_class(pt, "data.frame")
  expect_setequal(unique(pt$unit), seq_len(nrow(dec$quality)))
  expect_equal(nrow(pt), sum(dec$quality$n_discharges))
  expect_equal(pt$time_sec, (pt$discharge - 1) / sr)

  fs <- muFiringStats(dec)
  expect_equal(nrow(fs), nrow(dec$quality))
  expect_equal(fs$discharge_rate, dec$quality$discharge_rate)
  expect_true(all(fs$recruitment_sec <= fs$derecruitment_sec))
  expect_true(all(fs$cov_isi > 0 & fs$cov_isi < 0.5))     # regular firing
})

test_that("action-potential maps localise each unit on the grid", {
  maps <- muActionPotentialMaps(dec)
  expect_length(maps, nrow(dec$quality))
  m1 <- maps[[1]]
  expect_equal(dim(m1$map), sim$grid)
  expect_equal(ncol(m1$waveform), prod(sim$grid))
  expect_true(m1$peak_channel %in% seq_len(prod(sim$grid)))
  expect_equal(which.max(m1$amplitude), m1$peak_channel)
})

test_that("print methods report the objects", {
  expect_output(print(sim), "hdemg_sim")
  expect_output(print(dec), "motor unit")
})

# --- cheaper unit tests (small independent fixtures) --------------------------

test_that("make_hdemg_sim returns ground truth of the requested shape", {
  s <- make_hdemg_sim(n_units = 4, duration_sec = 2, grid = c(4, 5),
                      sampling_rate = 1000, seed = 3)
  expect_s3_class(s, "hdemg_sim")
  expect_equal(dim(s$signal), c(2000, 20))
  expect_length(s$spikes, 4)
  expect_length(s$muap, 4)
  expect_true(all(lengths(s$spikes) > 0))
  expect_true(all(s$discharge_rates > 3 & s$discharge_rates < 25))
  expect_true(all(vapply(s$muap, function(m) all(dim(m) == c(24, 20)),
                         logical(1))))
})

test_that("hdEMGDecompose validates its inputs", {
  expect_error(hdEMGDecompose(1:10), "matrix")
  expect_error(hdEMGDecompose("x", sampling_rate = 100), "matrix")
  expect_error(hdEMGDecompose(matrix(stats::rnorm(300), 100, 3)),
               "sampling_rate")
})

# --- regression tests for adversarial-review findings (WS5-10) ----------------

test_that("hdEMGDecompose rejects non-finite / non-scalar / non-positive knobs", {
  s <- make_hdemg_sim(n_units = 2, duration_sec = 2, grid = c(4, 4),
                      noise_sd = 0.005, seed = 1)
  # non-finite thresholds must fail fast, not corrupt the unit list
  expect_error(hdEMGDecompose(s, n_units = 2, extension_factor = 8,
                              pnr_threshold = NA), "pnr_threshold")
  expect_error(hdEMGDecompose(s, n_units = 2, extension_factor = 8,
                              sil_threshold = NaN), "sil_threshold")
  # non-finite / non-positive refractory_sec must not reach the index arithmetic
  expect_error(hdEMGDecompose(s, n_units = 2, extension_factor = 8,
                              refractory_sec = NA), "refractory_sec")
  expect_error(hdEMGDecompose(s, n_units = 2, extension_factor = 8,
                              refractory_sec = -0.02), "refractory_sec")
  # non-scalar knobs are rejected too
  expect_error(hdEMGDecompose(s, n_units = 2, extension_factor = 8,
                              pnr_threshold = c(20, 26)), "pnr_threshold")
  expect_error(hdEMGDecompose(s, n_units = 2, extension_factor = 8,
                              refractory_sec = c(0.02, 0.03)), "refractory_sec")
})

test_that("hdEMGDecompose rejects a non-finite signal and over-long extension", {
  s <- make_hdemg_sim(n_units = 2, duration_sec = 2, grid = c(4, 4), seed = 1)
  bad <- s$signal; bad[10, 3] <- NA
  expect_error(hdEMGDecompose(bad, sampling_rate = s$sampling_rate), "finite")
  expect_error(
    hdEMGDecompose(matrix(stats::rnorm(30), 10, 3), sampling_rate = 100,
                   extension_factor = 20),
    "extension_factor")
})

test_that("a bare matrix and an hdemg_sim decompose identically", {
  s <- make_hdemg_sim(n_units = 2, duration_sec = 3, grid = c(4, 4),
                      noise_sd = 0.004, seed = 5)
  d_mat <- hdEMGDecompose(s$signal, sampling_rate = s$sampling_rate,
                          n_units = 4, extension_factor = 10, seed = 1)
  d_sim <- hdEMGDecompose(s, n_units = 4, extension_factor = 10, seed = 1)
  expect_equal(nrow(d_mat$quality), nrow(d_sim$quality))
  if (length(d_mat$pulse_trains)) {
    expect_equal(d_mat$pulse_trains[[1]], d_sim$pulse_trains[[1]])
  }
})

test_that("a near-Gaussian recording yields no motor units", {
  set.seed(1)
  noise <- matrix(stats::rnorm(2000 * 9), 2000, 9)
  d <- hdEMGDecompose(noise, sampling_rate = 2048, n_units = 3,
                      extension_factor = 8)
  expect_s3_class(d, "hdemg_decomposition")
  expect_equal(nrow(muPulseTrains(d)), 0)
  expect_equal(nrow(muFiringStats(d)), 0)
  expect_length(muActionPotentialMaps(d), 0)
})

test_that("hdemgToPhysioExperiment round-trips through decomposition", {
  skip_if_not_installed("PhysioCore")
  skip_if_not_installed("S4Vectors")
  s <- make_hdemg_sim(n_units = 2, duration_sec = 3, grid = c(4, 4),
                      noise_sd = 0.004, seed = 1)
  pe <- hdemgToPhysioExperiment(s$signal, s$sampling_rate, grid = c(4, 4))
  expect_s4_class(pe, "PhysioExperiment")
  d <- hdEMGDecompose(pe, n_units = 4, extension_factor = 10, seed = 1)
  expect_equal(d$grid, c(4L, 4L))
  expect_error(
    hdemgToPhysioExperiment(s$signal, s$sampling_rate, grid = c(3, 3)),
    "rows")
})
