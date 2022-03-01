



#' @export
setGeneric('rmatrix', function(x, ...) standardGeneric('rmatrix'))




#' 
setMethod('rmatrix', 'numeric', 
		function(x, y=NULL, dist=runif, byrow = FALSE, dimnames = NULL, ...){
			
			x <- as.integer(x)
			# early exit if x has length 0
			if( length(x) == 0L )
				stop("KINOMO::rmatrix - invalid empty vector in argument `x`.")
			
			# check/ensure that 'dist' is a function.
			if( is.null(dist) ) dist <- runif
			if( isNumber(dist) ){
				os <- RNGseed()
				on.exit( RNGseed(os), add=TRUE)
				set.seed(dist)
				dist <- runif
			}
			if( !is.function(dist) )
				stop("KINOMO::rmatrix - invalid value for argument 'dist': must be a function [class(dist)='", class(dist), "'].")
			
			# if 'y' is not specified:
			if( is.null(y) ){
				
				if( length(x) == 1L ) y <- x # create a square matrix 
				else{ # assume x contains all dimensions (e.g. returned by dim())
					y <- x[2L]
					x <- x[1L]
				}
				
			}else{
				y <- as.integer(y)
				y <- y[1L] # only use first element
			}
			
			# build the random matrix using the distribution function
			matrix(dist(x*y, ...), x, y, byrow=byrow, dimnames=dimnames)	
		}
)



setMethod('rmatrix', 'ANY', 
		function(x, ...){
			rmatrix(x=dim(x), y=NULL, ...)
		}
)

