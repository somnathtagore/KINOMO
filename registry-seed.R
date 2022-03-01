

#' @include registry.R
#' @include KINOMOSeed-class.R
NULL

# create sub-registry for seeding methods
.registrySeed <- setPackageRegistry('seed', "KINOMOSeed"
		, description = "Initialization methods for KINOMO algorithms"
		, entrydesc = 'KINOMO seeding method')

KINOMOSeedInfo <- function(show=TRUE){
    obj <- .registrySeed
    if( show ) print(obj)
    invisible(obj)
}



#' 
KINOMOSeed <- function(name=NULL, ...){
	
	KINOMOGet('seed', name, ...)
	
}


getKINOMOSeed <- KINOMOSeed


existsKINOMOSeed <- function(name, exact=TRUE){	
	
	res <- !is.null( getKINOMOSeed(name, error=FALSE, exact=exact) )
	return(res)
	
}

# specific register method for registering KINOMOSeed objects
setMethod('KINOMORegister', signature(key='KINOMOSeed', method='missing'), 
		function(key, method, ...){
			KINOMORegister(name(key), key, ..., regname='seed')
		}
)


setKINOMOSeed <- function(..., overwrite=isLoadingNamespace(), verbose=TRUE){
	
	# wrap function method into a new KINOMOSeed object
	method <- KINOMOSeed(...)
	# register the newly created object
	res <- KINOMORegister(method, overwrite=overwrite, verbose=verbose)	
}

KINOMORegisterSeed <- setKINOMOSeed



removeKINOMOSeed <- function(name, ...){
	pkgreg_remove('seed', key=name, ...)
}

