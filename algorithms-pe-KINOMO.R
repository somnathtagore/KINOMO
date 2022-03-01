#' @include registry-algorithms.R
NULL



peKINOMO.objective <- function(fit, x, alpha, beta, ...)
{
	w <- .basis(fit)
	1/2 * sum( (x - fitted(fit))^2 )
		+ alpha * ( crossprod(w) - sum(w^2) )
		+ beta * sum(.coef(fit))
}

KINOMO_update.peKINOMO <- function(i, x, data, alpha, beta, ...){
	
	# retrieve each factor
	w <- .basis(data); h <- .coef(data);
	
	# At the first iteration initialise matrix M
	if( TRUE || i == 1 ){
		r <- ncol(w)
		M <- matrix(1, nrow=r, ncol=r) - diag(1, r)
		#staticVar('M', M, init=TRUE)
	}
	#else M <- staticVar('M')
	
	#precision threshold for numerical stability
	eps <- 10^-9
	
	# H_{i+1} = H_i ( W_i^T %*% V ) / ( W_i^T %*% W_i %*% H_i + beta)
	h <- h * crossprod(w, x) / ( crossprod(w) %*% h + beta)
	
	# W_{i+1} = W_i ( V %*% H_i^T ) / ( W_i %*% H_i %*% H_i^T + alpha W_i %*% M )
	w <- w * tcrossprod(x, h) / ( w %*% tcrossprod(h) + alpha * w %*% M )
	
	#return the modified data
	.basis(data) <- w; .coef(data) <- h;
	data
}

# register PE-KINOMO
KINOMOAlgorithm.peKINOMO <- setKINOMOMethod('pe-KINOMO', objective = peKINOMO.objective
		, model='KINOMOstd'
		, Update= KINOMO_update.peKINOMO
		, Stop='stationary')
