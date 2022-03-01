

#' @include rmatrix.R
NULL


KINOMOCheck <- function(method=NULL, rank=max(ncol(x)/5, 3), x=NULL, seed=1234, ...){
	
	# seed computation
	if( isNumber(seed) ){
		os <- RNGseed()
		on.exit( RNGseed(os), add=TRUE)
		set.seed(seed)
		seed <- NULL
	}
	if( is.null(x) ){
		x <- rmatrix(20, 10)
	}
	res <- KINOMO(x, rank, method, seed=seed, ...)
}
