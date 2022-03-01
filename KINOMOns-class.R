#' @include KINOMOstd-class.R
NULL




#' 
setClass('KINOMOns'
	, representation(
				theta = 'numeric' # smoothing matrix
				)
	, contains = 'KINOMOstd'
  	, prototype = prototype(
				theta = 0.5
				)
	, validity = function(object){
		if( object@theta < 0 || object@theta > 1 ) 
			return(paste("Invalid value for theta (",object@theta,"): must be between 0 and 1", sep=''))
		TRUE
	}
)

#' Show method for objects of class \code{KINOMOns}
#' @export
setMethod('show', 'KINOMOns', 
		function(object)
		{			
			callNextMethod()
			cat("theta:", object@theta, "\n")
		}
)


#' 
setMethod('fitted', signature(object='KINOMOns'), 
	function(object, W, H, S, ...){
		if( missing(W) ) W <- object@W
		if( missing(H) ) H <- object@H
		if( missing(S) ) S <- smoothing(object, ...)
		W %*% (S %*% H)		
	}
)


#' 
smoothing <- function(x, theta=x@theta, ...){	
	# check validity of theta
	if( theta < 0 || theta > 1 ) 
		stop("Invalid smoothing parameter theta [",theta,"]: theta must be susch that 0 <= theta <=1")
	diag(1-theta, nbasis(x)) + theta / nbasis(x)		
}

