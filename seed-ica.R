#' @include registry-seed.R
NULL


posICA <- function(object, x, ica.method=c('C', 'R'), ...){
			
	# perform ICA using the fastICA package
	if( !require.quiet('fastICA') )
		stop("Seeding method 'ica' requires package `fastICA` to be installed")
	requireNamespace('fastICA')
	ica.method <- match.arg(ica.method)
	res <- fastICA::fastICA(x, nbasis(object), method=ica.method, ...)
	
	# update the 'KINOMO' object
	.basis(object) <- pmax(res$S, .Machine$double.eps ); 
	.coef(object) <- pmax(res$A, .Machine$double.eps );	
	
	# return the updated object
	invisible(object)
	
} 

# Register positive ICA
setKINOMOSeed('ica', posICA, overwrite=TRUE)
