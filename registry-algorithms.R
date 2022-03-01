

#' @include registry.R
#' @include KINOMOStrategy-class.R
#' @include KINOMOStrategyFunction-class.R
#' @include KINOMOStrategyIterative-class.R
NULL

# create sub-registry for KINOMO algorithm
.registryAlgorithm <- setPackageRegistry('algorithm', "KINOMOStrategy"
		, description = "Algorithms to solve MF optimisation problems"
		, entrydesc = "KINOMO algorithm") 

KINOMOAlgorithmInfo <- function(show=TRUE){
    obj <- .registryAlgorithm
    if( show ) print(obj)
    invisible(obj)
}

# specific register method for registering KINOMOStrategy objects
setMethod('KINOMORegister', signature(key='KINOMOStrategy', method='missing'), 
		function(key, method, ...){
			KINOMORegister(name(key), key, ..., regname='algorithm')
		}
)


setKINOMOMethod <- function(name, method, ..., overwrite=isLoadingNamespace(), verbose=TRUE){
		
	# build call to KINOMOStrategy constructor
	call_const <- match.call(KINOMOStrategy)
	call_const[[1]] <- as.name('KINOMOStrategy')
	call_const$verbose <- NULL
	call_const$overwrite <- NULL
	# swap name and method if method is missing and name is a registered method
	if( missing(method) && !missing(name) && is.character(name) && existsKINOMOMethod(name) ){
		call_const$method <- name
		call_const$name <- NULL
	}
	# build the KINOMOStrategy object (in the parent frame to get the package slot right)
	e <- parent.frame()
	method <- eval(call_const, envir=e)
	# add to the algorithm registry
	res <- KINOMORegister(method, overwrite=overwrite, verbose=verbose)
	# return wrapper function invisibly
	wrap <- KINOMOWrapper(method)
}


KINOMORegisterAlgorithm <- setKINOMOMethod



NULL


setGeneric('canFit', function(x, y, ...) standardGeneric('canFit') )

#'    
setMethod('canFit', signature(x='KINOMOStrategy', y='character'),
		function(x, y, exact=FALSE){
			
			if( !exact ){
				
				# check for one model amongst all the models fittable by the strategy
				can <- if( length(mo <- modelname(x)) > 1 )
							sapply(mo, function(m) extends(y, m))
						else extends(y, mo)
				any(can)
				
			}else
				is.element(y, modelname(x))
		}
)
#' Tells if an KINOMO algorithm can fit the same class of models as \code{y}
setMethod('canFit', signature(x='KINOMOStrategy', y='KINOMO'),
		function(x, y, ...){
			canFit(x, modelname(y), ...)
		}
)
#' Tells if a registered KINOMO algorithm can fit a given KINOMO model
setMethod('canFit', signature(x='character', y='ANY'),
		function(x, y, ...){
			canFit(KINOMOAlgorithm(x), y, ...)
		}
)



selectKINOMOMethod <- function(name, model, load=FALSE, exact=FALSE, all=FALSE, quiet=FALSE){
	
	# lookup for an algorithm suitable for the given KINOMO model
	if( !isKINOMOclass(model) )
		stop("argument 'model' must be the name of a class that extends class 'KINOMO'")
	
	
	algo_list <- if( !missing(name) ){
				algo <- KINOMOAlgorithm(name)
				name(algo) 
			}else KINOMOAlgorithm()
	
	# lookup for all the algorithms that can fit the given model
	#NB: if only one model needs to be selected then first look for an exact fit as 
	# this would need to be done with exact=FALSE and TRUE anyways
	w <- sapply(algo_list, canFit, model, exact= if(all) exact else TRUE)	
	algo <- algo_list[w]
	
	# if no suitable algorithm was found, and an exact match is not required 
	# then look for other potential non-exact algorithms
	if( !all && !exact && length(algo) == 0 ){
		w <- sapply(algo_list, canFit, model, exact=FALSE)
		algo <- algo_list[w]
	}
	
	# return NULL if no algorithm was found
	if( length(algo) == 0L ){
		if( !quiet ) 
			stop("Could not find an KINOMO algorithm to fit model '", model, "'"
					, if( !missing(name) ) paste(" amongst ", str_out(algo_list, Inf)))
		return(NULL)
	}
	
	# if all=FALSE then try to choose the default algorithm if present in the list, or the first one
	res <- if( !all && length(algo) > 1L ){
				
				idx <- which( algo == KINOMO.getOption('default.algorithm') ) 
				if( !length(idx) ) idx <- 1L
				
				res <- algo[idx]
				if( !quiet ) 
					warning("Selected KINOMO algorithm '", res, "' amongst other possible algorithm(s): "
							, paste(paste("'", algo[-idx], "'", sep=''), collapse=", "))
				res
			}else # otherwise return all the algorithms
				algo
	
	# load the methods if required
	if( load ){
		if( length(res) > 1 ) sapply(res, KINOMOAlgorithm) else KINOMOAlgorithm(res)
	}
	else
		res	
}



getKINOMOMethod <- function(...) KINOMOGet('algorithm', ...)



#'  
KINOMOAlgorithm <- function(name=NULL, version=NULL, all=FALSE, ...){	
	
	# if one passes an KINOMOStrategy just returns it
	if( is(name, 'KINOMOStrategy') ) return(name)
	
	# force all=TRUE if type is provided
	if( !is.null(version) ) all <- TRUE
	
	# directly return the algorithm object if a key is supplied and all=FALSE
	if( !is.null(name) && !all ) return( getKINOMOMethod(name, ...) )
	
	# get all algorithms
	algo <- getKINOMOMethod(all=TRUE)
	# set names to match the primary key
	algo <- setNames(algo, sub("^\\.(.+#)?", '', algo))	
	# filter out hidden methods
	if( !all ) algo <- algo[!grepl("^\\.", algo)]
	# filter out methods not from the requested algorithm
	if( !is.null(name) ) algo <- algo[grepl(str_c("^", name), names(algo))]
	# filter out types
	if( !is.null(version)  ){
		type <- match.arg(version, c('R'))
		algo <- Filter( function(x) grepl(str_c("^\\.", version, '#'), x), algo)
	}
	
	# remove names if no arguments
	if( is.null(version) ) algo <- setNames(algo, NULL)
	# return the selected algorithm(s)
	algo
}



existsKINOMOMethod <- function(name, exact=TRUE){	
	
	!is.null( getKINOMOMethod(name, error=FALSE, exact=exact) )
	
}



removeKINOMOMethod <- function(name, ...){
	pkgreg_remove('algorithm', key=name, ...)
}




#' 
KINOMOWrapper <- function(method, ..., .FIXED=FALSE){
	
	# store original call
	.call <- match.call()
	
	# check that all arguments are named
	if( nargs() > 1L && any(names(.call)[-(1:2)]=='') )
		stop("Invalid call: all arguments must be named.")
	
	# store fixed arguments from default arguments
	.fixedargs <- 'method'
	.defaults <- names(.call)[-1L]
	.defaults <- .defaults[!.defaults %in% 'method']
	if( length(.defaults) ){
#		e <- parent.frame()
#		for(n in .defaults){
#			.call[[n]] <- eval(.call[[n]], envir=e)
#		}
		if( isTRUE(.FIXED) ) .fixedargs <- c(.fixedargs, .defaults)
		else if( is.character(.FIXED) ){
			.FIXED <- .FIXED[.FIXED %in% .defaults]
			.fixedargs <- c(.fixedargs, .FIXED)	
		}
	}
	# store in local environment
	.method <- method
	
	.checkArgs <- function(ca, args){
		# check for fixed arguments passed in the call that need
		# to be discarded
		nm <- names(ca)[-1L]
		if( any(fnm <- !is.na(pmatch(nm, .fixedargs))) ){
			warning("Discarding fixed arguments from wrapped call to ", .call[1L]
					, " [", str_out(nm[fnm], Inf), '].', immediate.=TRUE)
			ca <- ca[!c(FALSE, fnm)]
		}
		#
		
		# start with complete call
		.call <- ca
		# set values of wrapper default arguments if any
		if( length(.defaults) ){
			defaults <- args[.defaults]
			.call <- expand_list(ca, defaults, .exact=FALSE)
		}
		# change into a call to KINOMO
		.call[[1L]] <- as.name('KINOMO')
		.call[['method']] <- force(.method)
		as.call(.call)
	}
	
	# define wrapper function
	fwrap <- function(...){
		ca <- match.call()
		args <- formals()
		.call <- .checkArgs(ca, args)
		# eval in parent environment
		e <- parent.frame()
		eval(.call, envir=e)
	}
	
	# add default arguments to signature
	if( length(.defaults) ){
		formals(fwrap) <- expand_list(formals(fwrap), as.list(.call[.defaults]))
	}
	# add arguments from the KINOMO algorithm
	if( length(meth <- KINOMOFormals(.method)) ){
		formals(fwrap) <- expand_list(formals(fwrap), meth)
	}
	
	return( fwrap )
	
}

