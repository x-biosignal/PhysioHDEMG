test_that("readHDEMG round-trips an HDF5 HD-sEMG bundle into a PhysioExperiment", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("PhysioCore")
  skip_if_not_installed("S4Vectors")
  skip_if_not_installed("SummarizedExperiment")

  sim <- make_hdemg_sim(n_units = 3, duration_sec = 4, sampling_rate = 1024,
                        grid = c(4, 4), seed = 1)
  f <- tempfile(fileext = ".h5")
  on.exit(unlink(f), add = TRUE)
  rhdf5::h5createFile(f)
  rhdf5::h5write(sim$signal, f, "signal")            # n_time x n_channels
  rhdf5::h5write(sim$sampling_rate, f, "sampling_rate")
  rhdf5::h5write(as.integer(sim$grid), f, "grid")
  rhdf5::h5createGroup(f, "reference_mus")
  for (i in seq_along(sim$spikes)) {
    rhdf5::h5write(as.integer(sim$spikes[[i]]), f, paste0("reference_mus/mu_", i))
  }
  rhdf5::H5close()

  pe <- readHDEMG(f)
  expect_s4_class(pe, "PhysioExperiment")
  expect_equal(dim(SummarizedExperiment::assay(pe)), dim(sim$signal))
  expect_equal(PhysioCore::samplingRate(pe), sim$sampling_rate)
  refmu <- S4Vectors::metadata(pe)$reference_mupulses
  expect_length(refmu, length(sim$spikes))
  expect_equal(lengths(refmu), lengths(sim$spikes))
})

test_that("readHDEMG orients a channels x time signal to time x channels", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("PhysioCore")
  f <- tempfile(fileext = ".h5")
  on.exit(unlink(f), add = TRUE)
  sig_ct <- matrix(rnorm(6 * 400), nrow = 6, ncol = 400)   # 6 channels x 400 time
  rhdf5::h5createFile(f)
  rhdf5::h5write(sig_ct, f, "signal")
  rhdf5::h5write(500, f, "sampling_rate")
  rhdf5::H5close()
  pe <- readHDEMG(f, grid = c(2, 3))
  expect_equal(dim(SummarizedExperiment::assay(pe)), c(400L, 6L))  # -> time x channels
})

test_that("matchMotorUnits scores a decomposition against a reference (self-match = 1)", {
  sim <- make_hdemg_sim(n_units = 3, duration_sec = 6, sampling_rate = 1024, seed = 1)
  dec <- hdEMGDecompose(sim$signal, sampling_rate = sim$sampling_rate,
                        n_units = 5, pnr_threshold = 0, sil_threshold = 0)
  skip_if(length(dec$pulse_trains) < 1, "decomposition returned no units")

  m <- matchMotorUnits(dec, dec$pulse_trains, tol = 0L, max_lag = 0L)
  expect_s3_class(m, "mu_match")
  expect_true(all(c("ref_mu", "n_discharges", "matched_unit", "roa_lag_aware",
                    "roa_no_lag", "lag", "pnr_db", "sil") %in% names(m)))
  expect_equal(nrow(m), length(dec$pulse_trains))
  expect_true(all(m$roa_lag_aware == 1))             # each unit matches itself exactly
  expect_true(all(m$lag == 0L))
})

test_that("matchMotorUnits errors on an empty reference", {
  sim <- make_hdemg_sim(n_units = 2, duration_sec = 4, sampling_rate = 1024, seed = 2)
  dec <- hdEMGDecompose(sim$signal, sampling_rate = sim$sampling_rate,
                        n_units = 3, pnr_threshold = 0, sil_threshold = 0)
  expect_error(matchMotorUnits(dec, list()), "no reference")
})
