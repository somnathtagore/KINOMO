#' @include KINOMOstd-class.R
NULL



setClass('KINOMOOffset'
	, representation(
				offset = 'numeric' # offset vector
				)
	, contains = 'KINOMOstd'
  	, prototype=prototype(
  				offset = numeric()
				)
	
)

#' Show method for objects of class \code{KINOMOOffset}
#' @export
setMethod('show', 'KINOMOOffset', 
	function(object)
	{		
		callNextMethod()
		cat("offset: ")
		if( length(object@offset) > 0 ){
			cat('[', head(object@offset, 5)
				, if( length(object@offset) > 5 ) "..." else NULL
				, ']')
		}
		else cat('none')
		cat("\n")
	}
)


setMethod("initialize", 'KINOMOOffset', 
		function(.Object, ..., offset){			
			.Object <- callNextMethod()
			# correct the offset slot if possible
			if( missing(offset) ) offset <- numeric()
			if( !is.numeric(offset) ) stop("Unvalid value for parameter 'offset': a numeric vector is expected")			
			 
			# force length to be consistent with the factorization's dimension
			n <- nrow(.Object)
			if( n > 0 ) .Object@offset <- c( offset, rep(NA, max(0, n - length(offset))) )[1:n]
			
			# return the initialized valid object
			.Object
		}
)

#' @export
setGeneric('offset', package='stats')

#' 
setMethod('offset', signature(object='KINOMOOffset'), 
	function(object){
		object@offset
	}
)


setMethod('fitted', signature(object='KINOMOOffset'), 
	function(object, W, H, offset=object@offset){
		if( missing(W) ) W <- object@W
		if( missing(H) ) H <- object@H
		object@W %*% object@H + offset
	}
)


setMethod('rKINOMO', signature(x='KINOMOOffset', target='numeric'), 
function(x, target, ...){	
	
	# call the parent's 'rKINOMO' method to build a standard random KINOMO factorization
	res <- callNextMethod()
		
	#Vc# Initialize a random offset of length the number of genes
	res@offset <- runif(nrow(res), min=0, max=max(basis(res), coef(res)));
	
	# return the initialized KINOMOOffset object
	res
})
