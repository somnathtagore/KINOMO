#' @include KINOMOStrategy-class.R
NULL


#' 
setClass('KINOMOStrategyFunction'
	, representation(
		algorithm = 'function' # the function that implements the algorithm				
	)
	, contains = 'KINOMOStrategy'
)


setMethod('run', signature(object='KINOMOStrategyFunction', y='matrix', x='KINOMOfit'),
	function(object, y, x, ...){
		if( !is.function(fun <- algorithm(object)) )  
			stop("KINOMOStrategyFunction '", name(object), "': algorithm is not defined.")
		
		# run the function that defines the algorithm and return the result
		fun(y, x, ...)
	}
)

#' @export
KINOMOFormals.KINOMOStrategyFunction <- function(x, ...){
	args <- formals(x@algorithm)
	args[-(1:2)]
}
	

#' Returns the single R function that implements the KINOMO algorithm -- as stored in 
#' slot \code{algorithm}.
setMethod('algorithm', signature(object='KINOMOStrategyFunction'),
	function(object){
		slot(object, 'algorithm')
	}
)

#' Sets the function that implements the KINOMO algorithm, stored in slot \code{algorithm}. 
setReplaceMethod('algorithm', signature(object='KINOMOStrategyFunction', value='function'),
	function(object, value){
		slot(object, 'algorithm') <- value
		object
	}
)
