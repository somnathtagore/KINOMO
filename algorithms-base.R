

#' @include KINOMOstd-class.R
#' @include KINOMOOffset-class.R
#' @include KINOMOns-class.R
#' @include registry-algorithms.R
NULL




KINOMO_update.brunet_R <- function(i, v, x, eps=.Machine$double.eps, ...)
{
	# retrieve each factor
	w <- .basis(x); h <- .coef(x);
	
	# standard divergence-reducing KINOMO update for H
	h <- R_std.divergence.update.h(v, w, h)
	
	# standard divergence-reducing KINOMO update for W
	w <- R_std.divergence.update.w(v, w, h)
	
	#every 10 iterations: adjust small values to avoid underflow 
	if( i %% 10 == 0 ){
		#precision threshold for numerical stability
		#eps <- .Machine$double.eps
		h[h<eps] <- eps;
		w[w<eps] <- eps;
	}
	
	#return the modified model
	.basis(x) <- w; .coef(x) <- h;	
	return(x)
	
}


KINOMO_update.brunet <- function(i, v, x, copy=FALSE, eps=.Machine$double.eps, ...)
{
	# retrieve each factor
	w <- .basis(x); h <- .coef(x);
	# fixed terms
	nb <- nbterms(x); nc <- ncterms(x)
	
	# standard divergence-reducing KINOMO update for H	
	h <- std.divergence.update.h(v, w, h, nbterms=nb, ncterms=nc, copy=copy)
	
	# standard divergence-reducing KINOMO update for W
	w <- std.divergence.update.w(v, w, h, nbterms=nb, ncterms=nc, copy=copy)
	
	#every 10 iterations: adjust small values to avoid underflow
	# NB: one adjusts in place even when copy=TRUE, as 'h' and 'w' are local variables
	if( i %% 10 == 0 ){
		#eps <- .Machine$double.eps
		h <- pmax.inplace(h, eps, icterms(x))
		w <- pmax.inplace(w, eps, ibterms(x))
	}
	
	# update object if the updates duplicated the model
	if( copy ){		
		#return the modified model	
		.basis(x) <- w; 
		.coef(x) <- h;
	}
	return(x)
	
}



KINOMOAlgorithm.brunet_R <- setKINOMOMethod('.R#brunet'
		, objective='KL' 
		, Update=KINOMO_update.brunet_R
		, Stop='connectivity')


KINOMOAlgorithm.brunet <- setKINOMOMethod('brunet', '.R#brunet', Update=KINOMO_update.brunet)



KINOMOAlgorithm.KL <- setKINOMOMethod('KL'
		, objective='KL' 
		, Update=KINOMO_update.brunet
		, Stop='stationary')



KINOMO_update.lee_R <- function(i, v, x, rescale=TRUE, eps=10^-9, ...)
{
	# retrieve each factor
	w <- .basis(x); h <- .coef(x);	
	
	#precision threshold for numerical stability
	#eps <- 10^-9
	
	# compute the estimate WH
	#wh <- estimate(x)
	
	# euclidean-reducing KINOMO iterations	
	# H_au = H_au (W^T V)_au / (W^T W H)_au
	#h <- pmax(h * (t(w) %*% v),eps) / ((t(w) %*% w) %*% h + eps);
	h <- R_std.euclidean.update.h(v, w, h, eps=eps)
	
	# update H and recompute the estimate WH
	#metaprofiles(x) <- h
	#wh <- estimate(x)
	
	# W_ia = W_ia (V H^T)_ia / (W H H^T)_ia and columns are rescaled after each iteration	
	#w <- pmax(w * (v %*% t(h)), eps) / (w %*% (h %*% t(h)) + eps);
	w <- R_std.euclidean.update.w(v, w, h, eps=eps)
	#rescale columns TODO: effect of rescaling? the rescaling makes the update with offset fail
	if( rescale ) w <- sweep(w, 2L, colSums(w), "/", check.margin=FALSE)
	
	#return the modified model
	.basis(x) <- w; .coef(x) <- h;	
	return(x)
}



KINOMO_update.lee <- function(i, v, x, rescale=TRUE, copy=FALSE, eps=10^-9, weight=NULL, ...)
{
	# retrieve each factor
	w <- .basis(x); h <- .coef(x);
	# fixed terms
	nb <- nbterms(x); nc <- ncterms(x)
	
	#precision threshold for numerical stability
	#eps <- 10^-9
	
	# compute the estimate WH
	#wh <- estimate(x)
	
	# euclidean-reducing KINOMO iterations	
	# H_au = H_au (W^T V)_au / (W^T W H)_au
	h <- std.euclidean.update.h(v, w, h, eps=eps, nbterms=nb, ncterms=nc, copy=copy)
	# update original object if not modified in place
	if( copy ) .coef(x) <- h
	
	# W_ia = W_ia (V H^T)_ia / (W H H^T)_ia and columns are rescaled after each iteration	
	w <- std.euclidean.update.w(v, w, h, eps=eps, weight=weight, nbterms=nb, ncterms=nc, copy=copy)
	#rescale columns TODO: effect of rescaling? the rescaling makes the update with offset fail
	if( rescale ){
		w <- sweep(w, 2L, colSums(w), "/", check.margin=FALSE)
    }
	
	#return the modified model
	.basis(x) <- w; 	
	return(x)
}



KINOMOAlgorithm.lee_R <- setKINOMOMethod('.R#lee', objective='euclidean'
		, Update=KINOMO_update.lee_R
		, Stop='connectivity')	



KINOMOAlgorithm.lee <- setKINOMOMethod('lee', '.R#lee', Update=KINOMO_update.lee)

#

KINOMOAlgorithm.Frobenius <- setKINOMOMethod('Frobenius', objective='euclidean'
		, Update=KINOMO_update.lee
		, Stop='stationary')



KINOMO_update.euclidean_offset.h <- function(v, w, h, offset, eps=10^-9, copy=TRUE){
	.Call("offset_euclidean_update_H", v, w, h, offset, eps, copy, PACKAGE='KINOMO')
}
#' @export 
#' @rdname offset-KINOMO
KINOMO_update.euclidean_offset.w <- function(v, w, h, offset, eps=10^-9, copy=TRUE){
	.Call("offset_euclidean_update_W", v, w, h, offset, eps, copy, PACKAGE='KINOMO')
}
#' \code{KINOMO_update.offset_R} implements a complete single update step, 
#' using plain R updates.
#' @export 
#' @rdname offset-KINOMO
KINOMO_update.offset_R <- function(i, v, x, eps=10^-9, ...)
{	
	# retrieve each factor
	w <- .basis(x); h <- .coef(x);
	# retrieve offset and fill it if necessary (with mean of rows)
	off <- offset(x)
	if( i == 1 && length(off) == 0 )
		off <- rowMeans(v)
	
	#precision threshold for numerical stability
	#eps <- 10^-9
	
	# compute standard lee update (it will take the offset into account) without rescaling W's columns
	
	h <- R_std.euclidean.update.h(v, w, h, wh=w%*%h + off, eps=eps)
	w <- R_std.euclidean.update.w(v, w, h, wh=w%*%h + off, eps=eps)
	#x <- KINOMO_update.lee(i, v, x, rescale=FALSE, ...)
	
	# update the offset	
	# V0_i = V0_i ( sum_j V_ij ) / ( sum_j (V.off + W H)_ij )
	x@offset <- off * pmax(rowSums(v), eps) / (rowSums(w%*%h + off) + eps)
	
	#return the modified model
	.basis(x) <- w; .coef(x) <- h;
	return(x)
}
#' \code{KINOMO_update.offset} implements a complete single update step, 
#' using C++-optimised updates.
#' @export 
#' @rdname offset-KINOMO
KINOMO_update.offset <- function(i, v, x, copy=FALSE, eps=10^-9, ...)
{	
	# retrieve each factor
	w <- .basis(x); h <- .coef(x);
	# retrieve offset and fill it if necessary (with mean of rows)
	off <- offset(x)
	if( i == 1 && length(off) == 0 )
		off <- rowMeans(v)
	
	#precision threshold for numerical stability
	#eps <- 10^-9
	
	# compute standard offset updates
	h <- KINOMO_update.euclidean_offset.h(v, w, h, off, eps=eps, copy=copy)
	w <- KINOMO_update.euclidean_offset.w(v, w, h, off, eps=eps, copy=copy)
	
	# update the offset	
	# V0_i = V0_i ( sum_j V_ij ) / ( sum_j (V.off + W H)_ij )
	x@offset <- off * pmax(rowSums(v), eps) / (rowSums(w%*%h + off) + eps)	
	
	# update the original object if not modified in place
	if( copy ){ 
		.basis(x) <- w; 
		.coef(x) <- h;
	}
	return(x)
}



KINOMOAlgorithm.offset_R <- setKINOMOMethod('.R#offset', objective='euclidean'
		, model = 'KINOMOOffset'
		, Update=KINOMO_update.offset_R
		, Stop='connectivity')

# KINOMO with offset (optimised version)
#' @rdname offset-KINOMO
KINOMOAlgorithm.offset <- setKINOMOMethod('offset', '.R#offset', Update=KINOMO_update.offset)



#' @rdname nsKINOMO-KINOMO
KINOMO_update.ns <- function(i, v, x, copy=FALSE, ...)
{
	# retrieve and alter the factors for updating H
	S <- smoothing(x)
	w <- .basis(x)
	h <- .coef(x);
	
	# standard divergence-reducing update for H with modified W
	h <- std.divergence.update.h(v, w %*% S, h, copy=copy)
	
	# update H if not modified in place
	if( copy ) .coef(x) <- h
	
	# standard divergence-reducing update for W with modified H
	w <- std.divergence.update.w(v, w, S %*% h, copy=copy)
	
	# rescale columns of W
	w <- sweep(w, 2L, colSums(w), '/', check.margin=FALSE)
	
	#return the modified model
	.basis(x) <- w;
	return(x)
}
#' \code{KINOMO_update.ns_R} implements the same updates in \emph{plain R}.
#' 
#' @export
#' @rdname nsKINOMO-KINOMO
KINOMO_update.ns_R <- function(i, v, x, ...)
{
	# retrieve and alter the factors for updating H
	S <- smoothing(x)
	w <- .basis(x)
	#w <- metagenes(x) %*% smoothing(fit(x)); # W <- WS
	h <- .coef(x);
	
	# compute the estimate WH
	#wh <- estimate(x, W=w.init, H=h, S=S)
	
	# standard divergence-reducing update for H with modified W
	h <- R_std.divergence.update.h(v, w %*% S, h)
	
	# update H and recompute the estimate WH
	.coef(x) <- h
	# retrieve and alter the factors for updating W
	#w <- tmp;
	#h <- smoothing(fit(x)) %*% metaprofiles(x); # H <- SH
	#h <- S %*% h; # H <- SH
	
	# standard divergence-reducing update for W with modified H
	w <- R_std.divergence.update.w(v, w, S %*% h)
	
	# rescale columns of W
	w <- sweep(w, 2L, colSums(w), '/', check.margin=FALSE)
	
	#return the modified model
	.basis(x) <- w; #metaprofiles(x) <- h;
	return(x)
}



KINOMOAlgorithm.nsKINOMO_R <- setKINOMOMethod('.R#nsKINOMO', objective='KL'
		, model='KINOMOns'
		, Update=KINOMO_update.ns_R
		, Stop='connectivity')



KINOMOAlgorithm.nsKINOMO <- setKINOMOMethod('nsKINOMO', '.R#nsKINOMO', Update=KINOMO_update.ns)
