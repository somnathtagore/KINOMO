#' @include registry-seed.R
NULL



###% Seeding method for Nonnegative Matrix Factorization (NMF) algorithms.
###% 
###% @param object An instance of class \code{NMF} to seed
###% @param x The target matrix
###% @param method The method parameter passed to \code{fastICA}. Can either be 'R' or 'C' and
###% tells which implementation of fastICA to use (R code or C code).
###% @param ... extra parameters passed to \code{fastICA}
###% 
###% @return an updated version of \code{object}, where the matrix slots \code{W} and \code{H}
###% are set to the positive part of the IC of \code{x}.
###%  
posICA <- function(object, x, ica.method=c('C', 'R'), ...){
			
	# perform ICA using the fastICA package
	if( !require.quiet('fastICA') )
		stop("Seeding method 'ica' requires package `fastICA` to be installed")
    requireNamespace('fastICA')
	ica.method <- match.arg(ica.method)
	res <- fastICA::fastICA(x, nbasis(object), method=ica.method, ...)
	
	# update the 'NMF' object
	.basis(object) <- pmax(res$S, .Machine$double.eps ); 
	.coef(object) <- pmax(res$A, .Machine$double.eps );	
	
	# return the updated object
	invisible(object)
	
} 

# Register positive ICA
setNMFSeed('ica', posICA, overwrite=TRUE)
