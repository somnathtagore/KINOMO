

#' @include utils.R
NULL




#' 
syntheticKINOMO <- function(n, r, p, offset=NULL, noise=TRUE, factors=FALSE, seed=NULL){
	
	# set seed if necessary
	if( !is.null(seed) ){
		os <- RNGseed()
		on.exit( RNGseed(os) )
		set.seed(seed)
	}
	
	# internal parameters
	mu.W <- 1; sd.W <- 1
	if( isTRUE(noise) ){
		noise <- list(mean=0, sd=1)
	}else if( isNumber(noise) ){
		noise <- list(mean=0, sd=noise)
	}else if( is.list(noise) ){
		stopifnot( length(noise) == 2L )
		noise <- setNames(noise, c('mean', 'sd'))
	}else
		noise <- FALSE
	
	if( length(r) == 1 ){
		g <- rmultinom(1, p, rep(1, r))			
	}else{ # elements of r are the number of samples in each class 
		g <- r		
		p <- sum(r) # total number of samples
		r <- length(r) # number of class
	}
	
	# generate H
	H <- matrix(0, r, p)
	tmp <- 0
	for( i in 1:r ){
		H[i,(tmp+1):(tmp+g[i])] <- 1
		tmp <- tmp+g[i]
	} 	
	
	if( length(n) == 1 ){
		b <- rmultinom(1, n, rep(1, r))		
	}else{ # elements of n are the number of genes in each class 
		b <- n
		n <- sum(n)
	}
	
	# generate W
	W <- matrix(0, n, r)
	tmp <- 0
	for( i in 1:r ){		
		W[(tmp+1):(tmp+b[i]),i] <- abs(rnorm(b[i], mu.W, sd.W))
		tmp <- tmp + b[i]
	}	
	
	# build the composite matrix
	res <- W %*% H
	# add the offset if necessary
	if( !is.null(offset) ){
		if( length(offset) == 1L )
			offset <- rnorm(n, mean=0, sd=offset)
		
		stopifnot(length(offset)==n)
		res <- res +  offset
	}
	
	# add some noise if required
	if( !isFALSE(noise) )
		res <- pmax(res + rmatrix(res, dist=rnorm, mean=noise$mean, sd=noise$sd), 0)	
	
	# return the factors if required
	pData <- list(Group=factor(unlist(mapply(rep, 1:r, g, SIMPLIFY=FALSE))))
	fData <- list(Group=factor(unlist(mapply(rep, 1:r, b, SIMPLIFY=FALSE))))
	if( factors ) res <- list(res, W=W, H=H, offset=offset, pData=pData, fData=fData)
	
	# wrap results and expose relevant attributes
	ExposeAttribute(res, coefficients=H, basis=W, offset=offset
						, pData = pData, fData = fData
						, .VALUE=TRUE, .MODE='r')
	
}

