.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "PhysioHDEMG v", utils::packageVersion(pkgname),
    " - HD-sEMG motor-unit decomposition"
  )
}
