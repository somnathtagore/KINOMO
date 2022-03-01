

#' @include KINOMO-class.R
NULL


#' 
setGeneric('nneg', function(object, ...) standardGeneric('nneg'))




setMethod('nneg', 'matrix'
, function(object, method=c('pmax', 'posneg', 'absolute', 'min'), threshold=0, shift=TRUE){
	# match argument
	method <- match.arg(method)
	if( !is.numeric(threshold) || length(threshold) != 1L )
		stop("nneg - Invalid threshold value in argument `threshold` [",threshold,"]: must be a single numeric value.")
	if( threshold < 0 )
		stop("nneg - Invalid threshold value in argument `threshold` [",threshold,"]: must be nonnegative.")
	
	# 1. Transform if there is any negative entry
	m <- min(object)			
	if( m < 0 ){
		object <- 
		switch(method
		, pmax = pmax(object, 0)
		, posneg = rbind(pmax(object, 0), pmax(-object, 0))
		, absolute = pmax(abs(object), 0)
		, min = object - m
		, stop("KINOMO::nneg - Unexpected error: unimplemented transformation method '", method, "'.")
		)
	}

	if( threshold > 0 ){
		# 2. Apply threshold if any
		object <- pmax(object, threshold)
		
		# 3. Shifting: entries under threshold
		if( shift ) object[object<=threshold] <- 0
	}
	
	# return modified object
	object
}
)


setMethod('nneg', 'KINOMO', 
	function(object, ...){
		basis(object) <- nneg(basis(object), ...)
		object
	}
)


#' 
posneg <- function(...) nneg(..., method='posneg')


#' 
setGeneric('rposneg', function(object, ...) standardGeneric('rposneg'))


setMethod('rposneg', 'matrix'
, function(object, unstack=TRUE){
	
	# check that the number of rows is pair
	if( nrow(object) %% 2 != 0 )
		stop("rposneg - Invalid input matrix: must have a pair number of rows [",nrow(object),"].")
	n2 <- nrow(object)
	n <- n2/2
	if( unstack ) object <- object[1:n,,drop=FALSE] - object[(n+1):n2,,drop=FALSE]
	else object[(n+1):n2,] <- - object[(n+1):n2,,drop=FALSE]
	
	# return modified object
	object
}
)


setMethod('rposneg', 'KINOMO'
, function(object, ...){ 
	basis(object) <- rposneg(basis(object), ...)
	object
}
)



t.KINOMO <- function(x){
	# transpose and swap factors
	w <- t(basis(x))
	.basis(x) <- t(coef(x))
	.coef(x) <- w
	# return object
	x
}
