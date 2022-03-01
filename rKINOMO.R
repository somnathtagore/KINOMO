

#' @include KINOMOModel.R
NULL

.rKINOMO_fixed <- oneoffVariable('none')




#' 
setMethod('rKINOMO', signature(x='KINOMO', target='numeric'), 
	function(x, target, ncol=NULL, keep.names=TRUE, dist=runif){
		
		# store original dimnames
		if( keep.names ) dn <- dimnames(x)
		
		# valid parameter 'target'
		if( length(target) != 1 && length(target) != 2 )
			stop('KINOMO::rKINOMO - invalid target dimensions [length must be 1 or 2. Here length = ', length(target) ,']')
		if( any(is.na(target)) ) 
			stop('KINOMO::rKINOMO - invalid target dimensions [NA values in element(s): ', paste(which(is.na(target)), collapse=' and '), ']')		
		# shortcut for symetric case: provide only one dimension
		if( length(target) == 1L ){
			ncol <- if( !is.null(ncol) ){
				if( !is.numeric(ncol) || length(ncol) != 1 || is.na(ncol) )
					stop("KINOMO::rKINOMO - invalid argument `ncol`: must be a single numeric value")
				ncol
			}else target
			target <- c(target, ncol)
		}
		
		# retrieve dimension of the target matrix
		n <- target[1]; m <- target[2];
		# retrieve the factorization rank					
		r <- nbasis(x)
		
		## draw basis and coef matrices
		# interpret argument dist
		if( length(dist) == 0L ) dist <- runif
		if( is.character(dist) ){
			dist <- match.arg(dist, c('basis', 'coef'))
			dist <- setNames(list(runif), dist)
		}
		if( is.function(dist) ){
			dist <- list(basis = list(x=n, y=r, dist=dist)
					, coef = list(x=r, y=m, dist=dist))
		}else if( is.list(dist) ){
			if( !all(names(dist) %in% c('basis', 'coef')) ){
				dist <- list(basis=c(list(x=n, y=r), dist)
							, coef=c(list(x=r, y=m), dist))
			}else{
				if( !is.null(dist$basis) )
					dist$basis <- c(list(x=n, y=r), dist$basis)
				if( !is.null(dist$coef) )
					dist$coef <- c(list(x=r, y=m), dist$coef)
			}
		}
		
		fixed <- .rKINOMO_fixed()
		#Vc# Initialize random matrix: W
		# NB: this will keep the values of fixed basis terms
		if( !is.null(dist$basis) && !('basis' %in% fixed) ){
			basis(x) <- do.call('rmatrix', dist$basis);
		}
		#Vc# Initialize random matrix: H
		# NB: this will keep the values of fixed coef terms
		if( !is.null(dist$coef) && !('coef' %in% fixed) ){
			coef(x) <- do.call('rmatrix', dist$coef);
		}
		
		# if one needs to keep the names (possibly or reducing/increasing) 
		if( keep.names && !is.null(dn) )
			dimnames(x) <- list(dn[[1]][1:n], dn[[2]][1:m], dn[[3]][1:r])
		
		# return the modified object
		x
	}
)


#' 
setMethod('rKINOMO', signature(x='ANY', target='matrix'), 
	function(x, target, ..., dist=list(max=max(max(target, na.rm=TRUE), 1)), use.dimnames=TRUE){	
				
		# build a random KINOMO with the dimensions of the target matrix upper-bounded by the target's maximum entry.
		res <- rKINOMO(x, dim(target), ..., dist=dist)
		# compute the upper-bound of the random entries and enforce it if possible
		no.na <- abs(target[!is.na(target)])
		if( length(no.na) > 0 ){
			m <- max(no.na)
			basis(res) <- pmin(basis(res), m)
			coef(res) <- pmin(coef(res), m)
		}
		
		# set the dimnames from the target matrix if necessary
		if( use.dimnames )
			dimnames(res) <- dimnames(target)
		
		# return result
		res
	}
)
#' Shortcut for \code{rKINOMO(x, as.matrix(target))}.
setMethod('rKINOMO', signature(x='ANY', target='data.frame'),
	function(x, target, ...){
		rKINOMO(x, as.matrix(target), ...)
	}
)

#' 
setMethod('rKINOMO', signature(x='KINOMO', target='missing'), 
	function(x, target, ...){
		rKINOMO(x, c(nrow(x),ncol(x)), ...)
	}
)




#' 
setMethod('rKINOMO', signature(x='numeric', target='missing'),
	function(x, target, ..., W, H, dist=runif){
		
		# get fixed matrices to restore on exit:
		# one must enforce honouring the fixed matrices to prevent the call to 
		# rKINOMO from a sub-class method to change them.  
		of <- .rKINOMO_fixed()
		on.exit( .rKINOMO_fixed(of) )
		
		if( !missing(W) && missing(H) ){ # fixed basis matrix: x = n samples
			# one must not change the values H
			.rKINOMO_fixed('basis')
			x <- KINOMOModel(ncol(W), nrow(W), x, W=W, ...)
			dist <- list(coef=dist)
		}else if( missing(W) && !missing(H) ){ # fixed coef matrix: x = n features
			# one must not change the values H
			.rKINOMO_fixed('coef')
			x <- KINOMOModel(nrow(H), x, ncol(H), H=H, ...)
			dist <- list(basis=dist)
		}else if( !missing(W) && !missing(H) ){ # fixed basis and coef: x = rank
			# one must not change the values of W and H
			.rKINOMO_fixed(c('basis', 'coef'))
			x <- KINOMOModel(x, nrow(W), ncol(H), W=W, H=H, ...)
		}else
			stop("KINOMO::rKINOMO - Missing both arguments `W` and/or `H`: at least one of them must be specified.")
		
		rKINOMO(x, dist=dist) 
	}
)


setMethod('rKINOMO', signature(x='missing', target='missing'),
	function(x, target, ..., W, H){
		rKINOMO(min(ncol(W), nrow(H)), ..., W=W, H=H)
	}
)


#' 
setMethod('rKINOMO', signature(x='numeric', target='numeric'), 
		function(x, target, ncol=NULL, ..., dist=runif){		
			rKINOMO(KINOMOModel(x, target, ncol, ...), dist=dist)
		}
)
#' Generate a random formula-based KINOMO model, using the method 
#' \code{\link{KINOMOModel,formula,ANY-method}}.
setMethod('rKINOMO', signature(x='formula', target='ANY'), 
	function(x, target, ..., dist=runif){
		# missing target is NULL
		if( missing(target) ) target <- NULL
		rKINOMO(KINOMOModel(x, target, ...), dist=dist)
	}
)

