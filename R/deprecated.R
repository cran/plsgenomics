### deprecated.R  (2017-08)
###
###    Deprecated functions (that are now renamed) from the 'plsgenomics'
###    package
###
### Copyright 2017-08 Ghislain DURIF
###
###
### This file is part of the `plsgenomics' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

#' Deprecated function(s) in the 'plsgenomics' package
#' 
#' These functions are provided for compatibility with older version of
#' the 'plsgenomics' package.  They may eventually be completely
#' removed.
#' @rdname plsgenomics-deprecated
#' @name plsgenomics-deprecated
#' @param ... Parameters to be passed to the modern version of the function
#' @docType package
#' @export rirls.spls
#' @export rirls.spls.tune
#' @export rirls.spls.stab
#' @export m.rirls.spls
#' @export m.rirls.spls.tune
#' @export m.rirls.spls.stab
#' @export spls.adapt
#' @export spls.adapt.tune
#' @aliases
#' rirls.spls
#' rirls.spls.tune
#' rirls.spls.stab
#' m.rirls.spls
#' m.rirls.spls.tune
#' m.rirls.spls.stab
#' spls.adapt
#' spls.adapt.tune
#' @section Details:
#' \tabular{rl}{
#'     \code{rirls.spls} \tab is replaced by \code{\link{logit.spls}}\cr
#'     \code{rirls.spls.tune} \tab is replaced by \code{\link{logit.spls.cv}}\cr
#'     \code{rirls.spls.stab} \tab is replaced by \code{\link{logit.spls.stab}}\cr
#'     \code{m.rirls.spls} \tab is replaced by \code{\link{multinom.spls}}\cr
#'     \code{m.rirls.spls.tune} \tab is replaced by \code{\link{multinom.spls.cv}}\cr
#'     \code{m.rirls.spls.stab} \tab is replaced by \code{\link{multinom.spls.stab}}\cr
#'     \code{spls.adapt} \tab is replaced by \code{\link{spls}}\cr
#'     \code{spls.adapt.tune} \tab is replaced by \code{\link{spls.cv}}\cr
#' }
#'
m.rirls.spls <- function(...) {
     .Deprecated("multinom.spls", package="plsgenomics")
     multinom.spls(...)
}

m.rirls.spls.stab <- function(...) {
     .Deprecated("multinom.spls.stab", package="plsgenomics")
     multinom.spls.stab(...)
}

m.rirls.spls.tune <- function(...) {
     .Deprecated("multinom.spls.cv", package="plsgenomics")
     multinom.spls.cv(...)
}

rirls.spls <- function(...) {
     .Deprecated("logit.spls", package="plsgenomics")
     logit.spls(...)
}

rirls.spls.stab <- function(...) {
     .Deprecated("logit.spls.stab", package="plsgenomics")
     logit.spls.stab(...)
}

rirls.spls.tune <- function(...) {
     .Deprecated("logit.spls.cv", package="plsgenomics")
     logit.spls.cv(...)
}

spls.adapt <- function(...) {
     .Deprecated("spls", package="plsgenomics")
     spls(...)
}

spls.adapt.tune <- function(...) {
     .Deprecated("spls.cv", package="plsgenomics")
     spls.cv(...)
}
NULL
