

#' @include KINOMOstd-class.R
#' @include KINOMOns-class.R
#' @include KINOMOOffset-class.R
NULL


#' @export
#' @inline
setGeneric('KINOMOModel', function(rank, target=0L, ...) standardGeneric('KINOMOModel'))




#' 
setMethod('KINOMOModel', signature(rank='numeric', target='numeric'),
	function(rank, target, ncol=NULL, model='KINOMOstd', W, H, ..., force.dim=TRUE, order.basis=TRUE){
		
		if( is.null(model) ) model <- 'KINOMOstd'
		# check validity of the provided class
		if( !isClass(model) ) stop("KINOMOModel - Invalid model name: class '", model,"' is not defined.")
		if( !extends(model, 'KINOMO') ) stop("KINOMOModel - Invalid model name: class '", model,"' does not extend class 'KINOMO'.")
		
		# check the validity of the target
		if( length(target) == 0 ) stop('KINOMOModel - Invalid dimensions: `target` must be at least of length 1')
		if( length(target) > 2 ) stop('KINOMOModel - Invalid dimensions: `target` must be at most of length 2')
		if( !missing(ncol) && !is.null(ncol) && (!is.vector(ncol) || length(ncol) > 1 || !is.numeric(ncol) || ncol<0 ) )
			stop('KINOMOModel - Invalid dimensions: `ncol` must be a single nonnegative integer')
				
		# compute the target dimension
		target <- as.integer(target)
		n <- target[1]
		m <- if( length(target) == 2 ) target[2] 
			 else if( !missing(ncol) && !is.null(ncol) ) ncol
			 else if( !missing(H) ) ncol(H)
	 		 else n 
	 	if( n < 0 )
			stop("KINOMOModel - Invalid target number of rows: nonnegative value expected")
		if( m < 0 )
			stop("KINOMOModel - Invalid target number of columns: nonnegative value expected")
		# force rank to be an integer
		r <- as.integer(rank)
		
		# check the validity of the rank
		if( length(r) != 1 ) stop("Invalid argument 'rank': single numeric expected")
		if( r < 0 ) stop("KINOMOModel - Invalid argument 'rank': nonnegative value expected")
		
		# do not allow dimension incompatibility if required
		if( !force.dim && !missing(W) && !missing(H) && ncol(W) != nrow(H) ){
			stop('KINOMOModel - Invalid number of columns in the basis matrix [', ncol(W), ']: '
						, 'it should match the number of rows in the mixture coefficient matrix [', nrow(H), ']')
		}
				
		# build dummy compatible W and H if necessary
		W.was.missing <- FALSE
		if( missing(W) ){
			W <- matrix(as.numeric(NA), n, r)
			W.was.missing <- TRUE
		}
		else{
			if( is.vector(W) ) # convert numerical vectors into a matrix
				W <- matrix(W, n, r)
			else if( is.data.frame(W) ) # convert data.frame into matrix
				W <- as.matrix(W) 
			
			if( r == 0 ) r <- ncol(W)
			else if( r < ncol(W) ){
				if( !force.dim ){
					stop('KINOMOModel - Invalid number of columns in the basis matrix [', ncol(W), ']: ',
							'it should match the factorization rank [', r, ']')
				}
				warning("Objective rank is [",r,"] lower than the number of columns in W [",ncol(W),"]: "
						, "only the first ", r," columns of W will be used")
				W <- W[,1:r, drop=FALSE]				
			}
			else if( r > ncol(W) ){
				stop("KINOMOModel - Objective rank [",r,"] is greater than the number of columns in W [",ncol(W),"]")
			}
			
			# resolve consistency with target
			if( n == 0 ) n <- nrow(W)
			else if( n < nrow(W) ){
				if( !force.dim ){
					stop('KINOMOModel - Invalid number of rows in the basis matrix [', nrow(W), ']: '
							, 'it should match the target number of rows [', n, ']')
				}
				
				warning("KINOMOModel - Number of rows in target is lower than the number of rows in W [",nrow(W),"]: ",
							"only the first ", n," rows of W will be used")
				W <- W[1:n, , drop=FALSE]				
			}
			else if( n > nrow(W) ){
				stop("KINOMOModel - Number of rows in target [",n,"] is greater than the number of rows in W [",nrow(W),"]")
			}
		}
		
		if( missing(H) ) 
			H <- matrix(as.numeric(NA), ncol(W), m)
		else{
			# convert numerical vectors into a matrix
			if( is.vector(H) )
				H <- matrix(H, r, m)
			else if( is.data.frame(H) ) # convert data.frame into matrix
				H <- as.matrix(H)
			
			if( r == 0 ) r <- nrow(H)
			else if( r < nrow(H) ){
				if( !force.dim ){
					stop('KINOMOModel - Invalid number of rows in the mixture coefficient matrix [', nrow(H), ']: '
						, 'it should match the factorization rank [', r, ']')
				}
				warning("KINOMOModel - Objective rank [",r,"] is lower than the number of rows in H [",nrow(H),"]: "
								, "only the first ", r," rows of H  will be used")
				H <- H[1:r,, drop=FALSE]				
			}
			else if( r > nrow(H) ) stop("KINOMOModel - Objective rank [",r,"] is greater than the number of rows in H [",nrow(H),"]")
			# force dummy W to be at least compatible with H
			if( W.was.missing ) W <- matrix(as.numeric(NA), n, r)

			# resolve consistency with target
			if( m == 0 ) m <- ncol(H)
			else if( m < ncol(H) ){
				if( !force.dim ){
					stop('KINOMOModel - Invalid number of columns in the mixture coefficient matrix [', ncol(H), ']:'
						, ' it should match the target number of columns [', m, ']')
				}
				
				warning("KINOMOModel - Number of columns in target is lower than the number of columns in H [",ncol(H),"]:"
								, " only the first ", m," columns of H will be used")
				H <- H[, 1:m, drop=FALSE]				
			}
			else if( m > ncol(H) ){ 
				stop("KINOMOModel - Number of columns in target [",m,"]"
						," is greater than the number of columns in H [",ncol(H),"]")
			}
		}
		
		# check validity of matrices W and H (only if one of the target dimension is not null)
		if( n + m > 0 ){
			if( nrow(W) != n ) stop('KINOMOModel - Invalid number of rows for W: should match number of rows in target [', n, ']')
			if( ncol(W) != r ) stop('KINOMOModel - Invalid number of columns for W: should match factorization rank [', r, ']')
			if( nrow(H) != r ) stop('KINOMOModel - Invalid number of rows for H: should match factorization rank [', r, ']')
			if( ncol(H) != m ) stop('KINOMOModel - Invalid number of columns for H: should match number of columns in target [', m, ']')
		}
		
		# build and return a dummy KINOMO object
		KINOMO.debug('KINOMOModel', "Instantiate KINOMO model:", model)
		res <- new(model, ...)
		KINOMO.debug('KINOMOModel', "Set factors in model:", model)
		# set the dimnames if possible
		cW <- !is.null(colnames(W))
		rH <- !is.null(rownames(H))
		if( cW && !rH )# use colnames of W as basisnames
			rownames(H) <- colnames(W)
		else if( !cW && rH )# use rownames of H as basisnames
			colnames(W) <- rownames(H)
		else if( cW && rH ){# try to match names or use colnames of W (with a warning)
			
			# reorder as in the basis matrix if it makes sense, i.e. if the names are the same
			if( order.basis && !anyDuplicated(rownames(H)) && length(setdiff(rownames(H), colnames(W)))==0 ){ 
				H <- H[match(rownames(H), colnames(W)),]
			}
			else{
				rownames(H) <- colnames(W)
				warning("KINOMOModel - The rownames of the mixture matrix were set to match the colnames of the basis matrix")
			}
			
		}
		# set the basis and coef matrices
		.basis(res) <- W; .coef(res) <- H
		# check validity
		validObject(res)

		# return the model
		res
	}
)


#'  
setMethod('KINOMOModel', signature(rank='numeric', target='missing'),
		function(rank, target, ...){
			KINOMOModel(rank, 0L, ...)
		}
)


#' 
setMethod('KINOMOModel', signature(rank='missing', target='ANY'),
		function(rank, target, ...){
			KINOMOModel(0L, target, ...)
		}
)


setMethod('KINOMOModel', signature(rank='NULL', target='ANY'),
		function(rank, target, ...){
			KINOMOModel(0L, target, ...)
		}
)


#' 
setMethod('KINOMOModel', signature(rank='missing', target='missing'),
		function(rank, target, ...){
			# build an a priori empty model (extra args may provide the true dimension)
			# NB: do not allow dimension incompatibilities
			KINOMOModel(0L, 0L, ..., force.dim=FALSE)
		}
)


#' 
setMethod('KINOMOModel', signature(rank='numeric', target='matrix'),
		function(rank, target, ..., use.names=TRUE){
			# build an object compatible with the target's dimensions
			res <- KINOMOModel(rank, dim(target), ...)
			
			# try to set dimnames if it makes sense: 
			# set on target and not somehow already set on the result
			if( use.names && !is.null(dimnames(target)) ){
				dn <- dimnames(res)
				if( is.null(dn) )
					dn <- list(NULL, NULL, NULL)
				if( is.null(rownames(res)) && !is.null(rownames(target)) )
					dimnames(res) <- c(dimnames(target)[1], dn[2:3])
				if( is.null(colnames(res)) && !is.null(colnames(target)) )					
					dimnames(res) <- c(dimnames(res)[1], dimnames(target)[2], dimnames(res)[3])				
			}
			res
		}	
)



setMethod('KINOMOModel', signature(rank='matrix', target='matrix'),
		function(rank, target, ...){
			# use rank and target as W and H respectively
			# NB: do not allow dimension incompatibilities
			KINOMOModel(0L, 0L, W=rank, H=target, ..., force.dim=FALSE)
			
		}	
)



setMethod('KINOMOModel', signature(rank='data.frame', target='data.frame'),
	function(rank, target, ...){
		KINOMOModel(as.matrix(rank), as.matrix(target), ...)
	}
)



setMethod('KINOMOModel', signature(rank='matrix', target='ANY'),
		function(rank, target, ...){
			if( missing(target) ) target <- NULL
			# call KINOMOModel with swapping the arguments
			KINOMOModel(target, rank, ...)
			
		}	
)



parse_formula <- function(x){
	
	res <- list()
	# parse formula
	f <- as.character(x)
	hasResponse <- length(f) == 3L
	# response
	res$response <- hasResponse
	res$y <- if( hasResponse ) f[2L]
	# regressors
	reg <- if( hasResponse ) f[3L] else f[2L]
	res$x <- strsplit(reg, ' ')[[1]]
	res$n <- length(res$reg)
	# as a tring
	res$string <- paste(res$y, '~', reg, collapse='')
	
	res	
}

 the phenotypic data of x pData(x)
#' 
setMethod('KINOMOModel', signature(rank='formula', target='ANY'),
	function(rank, target, ..., data=NULL, no.attrib=FALSE){
		
		# missing target is NULL
		if( missing(target) ) target <- NULL
		
		# data is a model class name (passed from KINOMO)
		if( is.character(data) ){
			model <- data
			data <- NULL
		}else model <- NULL
		
		# parse formula
		f <- parse_formula(rank)
		enclos <- environment(rank)
		
		rank <- 0L
		if( is.vector(target) && is.numeric(target) ){
			rank <- target
			target <- NULL
		}
		
		# utility function to merge data and pData
		merge_pdata <- function(x, data){
			pd <- pData(x)
			if( length(pd) ){
				if( is.null(data) ) pd
				else{
					cbind(data, pd)
				}
			}else data
		}
		
		# determine formula data
		if( is.null(data) ){
			# target data.frame taken as data if a response variable if defined
			if( is.data.frame(target) && f$response ){
				data <- target
				target <- NULL
			}else if( is.environment(target) ){ # use target as enclosure
				enclos <- target
				target <- NULL 
			}
		}
		
		# determine target matrix:
		X <- 0L
		# if a response term is present, lookup target data in other arguments
		if( f$response ){
			X <- eval(parse(text=f$y), enclos)
			if( is.eset(target) && !identical(X, target) ){
				warning("Conflicting response term and target: the ExpressionSet in `target` will only be used for covariates.")
				data <- merge_pdata(target, data)
			}
		} 
		else if( is.null(target) ){
			# no response, no target: try ExpressionSet in data 
			if( is.eset(data) ){
				X <- exprs(data)
			}
		}else{
			X <- target
		}
		
		# merge data and pData from ExpressionSet target
		if( is.eset(X) ){
			data <- merge_pdata(X, data)
			X <- exprs(X)
		}
		
		r <- rank
		cterms <- bterms <- list()
		
		# dimensions are also inferred from the formula
		n <- if( identical(X, 0L) ) 0L else nrow(X)
		p <- if( identical(X, 0L) ) 0L else ncol(X)
		
		for( v in f$x ){
			if( grepl("^[0-9]+$", v) ){
				if( rank == 0L ){ # rank not specified in target 
					r <- as.numeric(v)
				}else{
					warning("KINOMO::KINOMOModel - Discarding rank specified in the formula [", v,"]:"
							, " using value specified in target rank instead [", rank, "].")
				}
			}else if( grepl("^[+-]$", v) ) next
			else {
				val <- eval(parse(text=v), data, enclos)
                .add_term <- function(v, val, type = NULL){
    				if( p==0L || length(val) ==  p || identical(type, 'coef') ){
    					cterms[[v]] <<- val
    					if( p==0L ) p <<- length(val)
    				}else if( n==0L || length(val) ==  n  || identical(type, 'basis') ){
    					bterms[[v]] <<- val
    					if( n==0L ) n <<- length(val)
    				}else
    					stop("Invalid", type," term '", v, "' length [", length(val), "]:"
                                , " length must either be the number of target columns [", p, "]"
    							, " or rows [", n, "]")
                }
                
                if( is.null(dim(val)) ) .add_term(v, val)
                else if( n == 0L || nrow(val) == n ){
                    lapply(1:ncol(val), function(i){
                        if( !is.null(cname <- colnames(val)[i]) && nzchar(cname) ) vname <- cname
                        else vname <- paste0(v, i)
                        .add_term(vname, val[, i], type = 'basis')   
                    })
                }else{
                    # special handling of data.frames: 
                    # -> coef terms are passed as column variables
                    if( is.data.frame(val) && (p == 0L || nrow(val) == p)){
                        val <- t(val)
                    } 
                    if( p == 0L || ncol(val) == p ){
                    lapply(1:nrow(val), function(i){
                        if( !is.null(cname <- rownames(val)[i]) && nzchar(cname) ) vname <- cname
                        else vname <- paste0(v, i)
                        .add_term(vname, val[i, ], type = 'coef')   
                    })
                    }else{
                        stop("Incompatible matrix-like term '", v, "' dimensions [", str_dim(val), "]:"
                                , " number of rows or columns must match the ones of the target matrix [", str_dim(X, dims = c(n, p)) ,"]")
                    }                
                }
			}
		}
		# try to fixup X if possible
		if( identical(X, 0L) ) X <- c(n, p)
		
		# call KINOMOModel with cterms
		if( hasArg(model) || is.null(model) ) object <- KINOMOModel(r, X, ...)
		else object <- KINOMOModel(r, X, ..., model=model)
		# set fixed basis terms
		if( length(bterms) ){
			bterms(object) <- as.data.frame(bterms)	
		}
		# set fixed coef terms
		if( length(cterms) ){
			cterms(object) <- as.data.frame(cterms)
		}
		
		# valid object
		validObject(object)
		# attach formula data
		if( !no.attrib ){
			attr(object, 'target') <- X
			attr(object, 'formula') <- f
		}
		# return object
		object
	}
)



#' 
KINOMOModels <- function(builtin.only=FALSE){
	
	if( builtin.only ) return( .KINOMO.Models.Builtin )
	
	# return all subclasses of class 'KINOMO' (minus class 'KINOMOfit' and its subclasses)
	models <- names(methods::getClass('KINOMO')@subclasses)
	models.wraps <- c('KINOMOfit', names(methods::getClass('KINOMOfit')@subclasses))
	return( models[!is.element(models, models.wraps)] )
	
}

###% Initialization function for KINOMO models
.KINOMO.Models.Builtin <- NULL
.init.KINOMO.models <- function(){	
	.KINOMO.Models.Builtin <<- KINOMOModels()
}

