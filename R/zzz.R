.onAttach <- function(libname, pkgname) {
     packageStartupMessage("For any news related to the 'plsgenomics' package (update, corrected bugs), please check http://thoth.inrialpes.fr/people/gdurif/")
     packageStartupMessage("C++ based sparse PLS routines will soon be available on the CRAN in the new 'fastPLS' package.")
}

.onLoad <- function(libname, pkgname) {
     ## setting BLAS and OPENMP number of threads (modified in each function)
     blas_set_num_threads(1)
     omp_set_num_threads(1)
}
