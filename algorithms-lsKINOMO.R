

#' @include registry-algorithms.R
NULL


KINOMO_update.lsKINOMO <- function(i, X, object, weight, eps=10^-9, ...)
{
	if( i == 1 ){# pre-compute weighted target matrix
		staticVar('wX', X * weight, init=TRUE)
	}
	
	# retrieve weighted target matrix
	wX <- staticVar('wX')
	
	# retrieve each factor
	w <- .basis(object); h <- .coef(object);	
	
	# compute the estimate WH
	wh <- fitted(object) * weight
	
	# euclidean-reducing KINOMO iterations	
	# H_au = H_au (W^T V/sigma)_au / (W^T (W H)/sigma)_au
	h <- KINOMO_update.euclidean.h_R(wX, w, h, wh=wh, eps=eps)	
	
	# update H and recompute the estimate WH
	.coef(object) <- h;
	wh <- fitted(object) * weight
	
	# W_ia = W_ia (V/sigma H^T)_ia / ((W H)/sigma H^T)_ia and columns are rescaled after each iteration	
	w <- KINOMO_update.euclidean.w_R(wX, w, h, wh=wh, eps=eps)	
	
	#return the modified data
	.basis(object) <- w	
	return(object)
}

#' \code{wrss} implements the objective function used by the LS-KINOMO algorithm.
#' 
#' @rdname lsKINOMO-KINOMO
wrss <- function(object, X, weight){
	sum( ((X - fitted(object)) * weight)^2 )/2
}


KINOMOAlgorithm.lsKINOMO <- setKINOMOMethod('ls-KINOMO', objective=wrss
			, Update=KINOMO_update.lsKINOMO
			, Stop='stationary')
	
# Unit test for the LS-KINOMO algorithm
runit.lsKINOMO <- function(){
	
	requireNamespace('RUnit')
	set.seed(12345)
	X <- rmatrix(100,20)
	
	res <- KINOMO(X, 3, 'ls-KINOMO', weight=1, seed=1)	
	res2 <- KINOMO(X, 3, '.R#lee', rescale=FALSE, seed=1, .stop=KINOMO.stop.stationary)
	tol <- 10^-14
	RUnit::checkTrue( KINOMO.equal(res, res2, identical=FALSE, tol=tol ), paste("LS-KINOMO with weight = 1 and .R#Lee (no scale + stationary) give identical results at tolerance=", tol))	
	
}
