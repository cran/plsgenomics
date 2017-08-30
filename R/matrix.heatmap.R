### Copyright 2016-06 Ghislain DURIF
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

#' @title
#' Heatmap visualization for matrix
#' @aliases matrix.heatmap
#' 
#' @description 
#' Visualization of matrix entries in heatmap format, the color scale 
#' depends on the numerical values.
#' 
#' @details
#' The function \code{matrix.heatmap} is a wrapper for the function 
#' \code{\link[fields]{image.plot}} from the 'fields' package.
#' 
#' @param mat the matrix to visualize
#' @param ... any argument that could be pass to the functions 
#' \code{\link[fields]{image.plot}} or \code{\link[graphics]{image}}.
#' 
#' @return No return, just plot the heatmap in the current graphic window.
#' 
#' @author
#' Ghislain Durif (\url{http://thoth.inrialpes.fr/people/gdurif/}).
#' 
#' @seealso \code{\link{logit.spls}}, \code{\link{stability.selection}}, 
#' \code{\link{stability.selection.heatmap}}
#' 
#' @examples
#' ### load plsgenomics library
#' library(plsgenomics)
#' 
#' ### generate a matrix
#' A = matrix(runif(10*10), ncol=10)
#' 
#' ### heatmap of estimated probabilities
#' matrix.heatmap(A)
#' 
#' @export
matrix.heatmap <- function(mat, ...) {
     image.plot(t(apply(mat, 2, rev)), xaxt="n", yaxt="n", ...)
}
