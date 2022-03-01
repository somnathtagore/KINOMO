# Algorithm for Nonnegative Matrix Factorization: Local KINOMO (LKINOMO)
#
# @author Renaud Gaujoux
# @created 21 Jul 2009

#' @include registry-algorithms.R
NULL



KINOMO_update_R.lKINOMO <- function(i, v, data, ...){
	
	# retrieve each factor
	w <- .basis(data); h <- .coef(data);
	
	# update H 
	h <- sqrt( h * crossprod(w, v / (w %*% h)) )
	
	# update W using the standard divergence based update
	w <- R_std.divergence.update.w(v, w, h, w %*% h)
	
	# scale columns of W
	w <- sweep(w, 2L, colSums(w), "/", check.margin=FALSE)	
	
	#every 10 iterations: adjust small values to avoid underflow 
	if( i %% 10 == 0 ){
		#precision threshold for numerical stability
		eps <- .Machine$double.eps
		h[h<eps] <- eps;
		w[w<eps] <- eps;
	}
		
	# return updated data	
	.basis(data) <- w; .coef(data) <- h
	return(data)
}

KINOMO_update.lKINOMO <- function(i, v, data, ...){
	
	# retrieve each factor
	w <- .basis(data); h <- .coef(data);
	
	# update H 
	h <- sqrt( h * crossprod(w, v / (w %*% h)) )
	
	# update W using the standard divergence based update
	w <- std.divergence.update.w(v, w, h)
	
	# scale columns of W
	w <- apply(w, 2, function(x) x/sum(x))
	
	#every 10 iterations: adjust small values to avoid underflow 
	if( i %% 10 == 0 ){
		#precision threshold for numerical stability
		eps <- .Machine$double.eps
		h[h<eps] <- eps;
		w[w<eps] <- eps;
	}
	
	# return updated data	
	.basis(data) <- w; .coef(data) <- h
	return(data)
}


