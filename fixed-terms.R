

#' @include KINOMO-class.R
#' @include KINOMOstd-class.R
NULL


setMethod('c', 'KINOMO',
	function(x, ..., margin=3L, recursive=FALSE){
		
		y <- ..1
		if( is.matrix(y) ){
			
			if( missing(margin) ){
				if( nrow(y) == nrow(x) ){
					if( ncol(y) == ncol(x) ){
						warning("KINOMO::`c` - Right argument match both target dimensions: concatenating basis columns."
								, " Use `margin=4L` to concatenate coefficient rows.")
					}
					margin <- 3L
				}else if( ncol(y) == ncol(x) ){
					margin <- 4L
				}else{
					stop("KINOMO::`c` - Incompatible argument dimensions: could not infer concatenation margin.")
				}
			}
			
			if( margin == 1L ){ # extend basis vectors
				
				if( nbterms(x) ){ # cannot extend models with fixed basis terms
					stop("KINOMO::`c` - Could not extend basis vectors:"
							, " KINOMO model has fixed basis terms [", nbterms(x), "]")
				}
				if( ncol(y) != nbasis(x) ){
					stop("KINOMO::`c` - Could not extend basis vectors:"
							, " incompatible number of columns [", nbasis(x), '!=', ncol(y), "].")
				}
				
				# extend basis vectors
				basis(x) <- rbind(basis(x), y)
				
			} else if( margin == 2L ){ # extend basis profiles
				
				if( ncterms(x) ){ # cannot extend models with fixed coef terms
					stop("KINOMO::`c` - Could not extend basis profiles:"
							, " KINOMO model has fixed coefficient terms [", ncterms(x), "]")
				}
				if( nrow(y) != nbasis(x) ){
					stop("KINOMO::`c` - Could not extend basis profiles:"
							, " incompatible number of rows [", nbasis(x), '!=', nrow(y), "].")
				}
				
				# extend basis profiles
				coef(x) <- cbind(coef(x), y)
				
			} else if( margin == 3L ){ # add basis vectors
				
				if( nrow(y) != nrow(x) ){
					stop("KINOMO::`c` - Could not concatenate basis vectors:"
						, " incompatible number of rows [", nrow(x), '!=', nrow(y), "].")
				}
				
				# bind basis terms
				.basis(x) <- cbind(basis(x), y)
				dn <- colnames(.basis(x))
				# bind dummy coef
				.coef(x) <- rbind(coef(x), matrix(NA, ncol(y), ncol(x)))
				basisnames(x) <- dn
				
			} else if( margin == 4L ){ # add basis profiles
				
				if( ncol(y) != ncol(x) ){
					stop("KINOMO::`c` - Could not concatenate basis profiles:"
						, " incompatible number of columns [", ncol(x), '!=', ncol(y), "].")
				}
				
				# bind coef terms
				.coef(x) <- rbind(coef(x), y)
				dn <- rownames(.coef(x))
				# bind dummy basis
				.basis(x) <- cbind(basis(x), matrix(NA, nrow(x), nrow(y)))
				basisnames(x) <- dn
				
			}else{
				stop("KINOMO::`c` - Invalid concatenation margin: should be either"
								, " 1L (basis rows), 2L (coef columns), 3L (basis vectors/columns) or 4L (basis profiles/coef rows).")
			}
		}else if( is.KINOMO(y) ){
			# check dimensions
			if( nrow(x) != nrow(y) )
				stop("KINOMO::`c` - Could not concatenate KINOMO objects:"
					, " incompatible number of rows [", nrow(x), '!=', nrow(y), "]")
			if( ncol(x) != ncol(y) )
				stop("KINOMO::`c` - Could not concatenate KINOMO objects:"
					, " incompatible number of columns [", ncol(x), '!=', ncol(y), "]")
			.basis(x) <- cbind(basis(x), basis(y))
			.coef(x) <- rbind(coef(x), coef(y))
		}else{
			stop("KINOMO::`c` - Concatenation of an KINOMO object with objects of class '", class(y), "' is not supported.")
		}
		
		# return augmented object
		x
	}
)


fterms <- function(value){

	res <- list(n=0L, terms=NULL, df=NULL, i=integer())
	if( !length(value) ) return(res)
	
	# convert into a data.frame
	if( is.factor(value) ) value <- data.frame(Group=value)
	else if( is.numeric(value) ) value <- data.frame(Var=value)
	else if( !is.data.frame(value) ) value <- as.data.frame(value)
	
	res$n <- length(value)
    res$df <- value
	# generate fixed term matrix
	terms <- model.matrix(~ -1 + ., data=value)
	res$terms <- terms
	# build indexes
	res$i <- 1:ncol(terms)
		
	res
}


setGeneric('bterms<-', function(object, value) standardGeneric('bterms<-'))

setReplaceMethod('bterms', signature('KINOMOstd', 'ANY'),
	function(object, value){
		
		if( nterms(object) ){
			stop("Cannot set fixed basis terms on an object that already has fixed terms:",
					" these can be set only once and before setting any fixed coefficient term",
					" [coef=", ncterms(object), ", basis=", nbterms(object), "].")
		}
		# build terms
		t <- fterms(value)
		
		if( !t$n ) return(object)
		
		# check dimension
		if( nrow(t$terms) != nrow(object) ){
			stop("Invalid fixed basis terms: all terms should have length the number of target rows"
					, "[terms=", nrow(t$terms), " != ", nrow(object), "=target]")
		}
		
		# set data
		object@bterms <- t$df
		# set indexes
		i <- t$i
		nv <- nbasis(object)
		object@ibterms <- nv + i
		# set terms
		object <- c(object, t$terms, margin=3L)
		
		object
	}
)

setGeneric('cterms<-', function(object, value) standardGeneric('cterms<-'))

setReplaceMethod('cterms', signature('KINOMOstd', 'ANY'),
	function(object, value){
		
		if( ncterms(object) ){
			stop("Cannot set fixed coef terms on an object that already has fixed coef terms:",
				" these can be set only once", 
				" [coef=", ncterms(object), ", basis=", nbterms(object), "].")
		}
		# build terms
		t <- fterms(value)
		
		if( !t$n ) return(object)
		
		# check dimension
		if( nrow(t$terms) != ncol(object) ){
			stop("Invalid fixed coefficient terms: all terms should have length the number of target columns"
					, "[terms=", nrow(t$terms), " != ", ncol(object), "=target]")
		}
		
		# transpose term matrix
		t$terms <- t(t$terms)
		# set data
		object@cterms <- t$df
		# set indexes
		i <- t$i
		nv <- nbasis(object)
		object@icterms <- nv + i
		# set terms
		object <- c(object, t$terms, margin=4L)
					
		object
	}
)


setGeneric('ibterms', function(object, ...) standardGeneric('ibterms') )

setMethod('ibterms', 'KINOMO', 
	function(object, ...){
		stop("KINOMO::ibterms is a pure virtual method of interface 'KINOMO'."
			," It should be overloaded in class '", class(object),"'.")
	}
)

setMethod('ibterms', 'KINOMOstd', 
	function(object){
		object@ibterms
	}
)

setGeneric('icterms', function(object, ...) standardGeneric('icterms') )

setMethod('icterms', 'KINOMO', 
	function(object, ...){
		stop("KINOMO::icterms is a pure virtual method of interface 'KINOMO'."
				," It should be overloaded in class '", class(object),"'.")
	}
)

setMethod('icterms', 'KINOMOstd', 
	function(object){
		object@icterms
	}
)

iterms <- function(object, ...){
	c(ibterms(object), icterms(object))
}


nterms <- function(object){
	length(ibterms(object)) + length(icterms(object))
}

nbterms <- function(object){
	length(ibterms(object))
}

ncterms <- function(object){
	length(icterms(object))
}


bterms <- function(object){
	object@bterms		
}

cterms <- function(object){
	object@cterms	
}


ibasis <- function(object, ...){
	i <- 1:nbasis(object)
	if( length(idx <- ibterms(object, ...)) ) i[-idx]
	else i
}

icoef <- function(object, ...){
	i <- 1:nbasis(object)
	if( length(idx <- icterms(object, ...)) ) i[-idx]
	else i
}

#' @export
t.KINOMOstd <- function(x){
    # transpose and swap factors
    x <- t.KINOMO(x)
    # swap fixed terms
    bt <- bterms(x)
    ibt <- ibterms(x)
    x@bterms <- cterms(x)
    x@ibterms <- icterms(x)
    x@cterms <- bt
    x@icterms <- ibt
    
    # returns
    x
}
