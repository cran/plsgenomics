### utils.R  (2017-08)
###
###    Utility functions
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


safeExp <- function(x) {
     if(x < -27) {
          return(0)
     } else if(x > 27) {
          return(1)
     } else {
          return(exp(x))
     }
}

safeExpMat <- function(X) {
     return(apply(X, c(1,2), safeExp))
}

safeSum <- function(x) {
     if(sum(x) == 0) {
          return(1E-7)
     } else {
          return(sum(x))
     }
}

softMax <- function(eta) {
     t(apply(safeExpMat(eta), 1, function(x) x/safeSum(x)))
}
