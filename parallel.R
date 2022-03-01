

#' @include utils.R
#' @import foreach
#' @import doParallel
NULL


getMaxCores <- function(limit=TRUE){
	#ceiling(parallel::detectCores()/2)
	nt <- n <- parallel::detectCores()
	# limit to number of cores specified in options if asked for
	if( limit ){
		if( !is.null(nc <- getOption('cores')) ) n <- nc # global option
		else if( !is.null(nc <- KINOMO.getOption('cores')) ) n <- nc # KINOMO-specific option
		else if( n > 2 ) n <- n - 1L # leave one core free if possible
	}
	# forces limiting maximum number of cores to 2 during CRAN checks
	if( n > 2 && isCHECK() ){
		message("# NOTE - CRAN check detected: limiting maximum number of cores [2/", nt, "]")
		n <- 2L
	}
	n
}


registerDoBackend <- function(object, ...){

	# restore old backend data in case of an error
	old <- getDoBackend()
	on.exit( setDoBackend(old) )
	
	# get old foreach backend object
	ob <- ForeachBackend()
	
	# register new backend: call the register method
	b <- ForeachBackend(object, ...)
	res <- register(b)
	
	# cancel backend restoration
	on.exit()
	# call old backend cleanup method
	doBackendCleanup(ob)
	
	# return old backend
	invisible(ob)
}


getDoBackend <- function(){
    fe_ns <- asNamespace('foreach')
	fe <- ns_get('.foreachGlobals', fe_ns)
	if( !exists("fun", where = fe, inherits = FALSE) )
		return(NULL)
	
    getDoPar <- ns_get('getDoPar', fe_ns)
	c(getDoPar() # this returns the registered %dopar% function + associated data
		# -> add info function from foreach internal environment
		, info= if( exists("info", where = fe, inherits = FALSE) ){
					get('info', fe, inherits=FALSE) 
				}else{
					function(data, item) NULL
				}
		, cleanup = if( exists("cleanup", where = fe, inherits = FALSE) ){
			get('cleanup', fe, inherits=FALSE)
		}
	)
}

getDoBackendInfo <- function(x, item){
    if( is.function(x$info) ) x$info(x$data, item)
}
getDoBackendName <- function(x){
    getDoBackendInfo(x, 'name')
}


setDoBackend <- function(data, cleanup=FALSE){
	
	# get old backend data
	ob <- getDoBackend()
	ofb <- ForeachBackend()
	# cleanup old backend if requested
	if( cleanup ){
		doBackendCleanup(ofb)
	}
	
	if( !is.null(data) ){
		bdata <- data
		if( is.backend(data) ) data <- data[!names(data) %in% c('name', 'cleanup')]
		do.call('setDoPar', data)
		setBackendCleanup(bdata)
	}else{
		do.call('setDoPar', list(NULL))
		fe <- ns_get('.foreachGlobals', 'foreach')
		if (exists("fun", envir = fe, inherits = FALSE))
			remove("fun", envir = fe)
		setBackendCleanup(NULL)
	}
	# return old backend
	invisible(ob)
}

# setup cleanup procedure for the current backend
setBackendCleanup <- function(object, fun, verbose=FALSE){
	
	fe <- ns_get('.foreachGlobals', 'foreach')
	name <- getDoParName()
	if( !is.null(fun <- object$cleanup) ){
		if( verbose ) message("# Registering cleaning up function for '", name, "'... ", appendLF=FALSE)
		assign('cleanup', fun, fe)
		if( verbose ) message("OK")
	}else if (exists("cleanup", envir = fe, inherits = FALSE)){
		if( verbose ) message("# Removing cleaning up function for '", name, "'... ", appendLF=FALSE)
		remove("cleanup", envir = fe)
		if( verbose ) message("OK")
	}
	invisible(object)
}

# run cleanup procedure for a given backend object
doBackendCleanup <- function(object, ..., run=TRUE, verbose=FALSE){
	
	name <- object$name
	if( !is.null(fun <- object$cleanup) ){
		if( verbose ) message("# Cleaning up '", name, "'... ", appendLF=FALSE)
		res <- try(fun(), silent=TRUE) 
		if( verbose ) message(if( is(res, 'try-error') ) 'ERROR' else 'OK')
		if( isTRUE(res) ) object$cleanup <- NULL
		if( verbose ) message('OK', if( !is.null(res) ) str_c(' [', res,']'))
	}
	invisible(object)
}


register <- function(x, ...){
	UseMethod('register', x)
}
#' @export
register.foreach_backend <- function(x, ...){
	
	be <- x$name
	# For everything except doSEQ:
	# require definition package (it is safer to re-check)
	if( be != 'doSEQ' ){
		if( !require.quiet(be, character.only=TRUE) )
			stop("Package '", be, "' is required to use foreach backend '", be, "'")
	}
	
	regfun <- .foreach_regfun(x$name)
	res <- 
	if( length(formals(regfun)) > 0L ) do.call(regfun, c(x$data, ...))
	else regfun()
	# throw an error if not successful (foreach::setDoPar do not throw errors!!)
	if( is(res, 'simpleError') ) stop(res)
	# set cleanup procedure if any
	setBackendCleanup(x)
	# return result
	invisible(res)
}


#' @rdname foreach
setGeneric('ForeachBackend', function(object, ...) standardGeneric('ForeachBackend'))

#' Default method defined to throw an informative error message, when no other
#' method was found.
setMethod('ForeachBackend', 'ANY', 
	function(object, ...){
		if( is.backend(object) ){
			# update arg list if necessary
			if( nargs() > 1L )	object$data <- list(...)
			object
		}else if( is(object, 'cluster') )
			selectMethod('ForeachBackend', 'cluster')(object, ...)
		else
			stop("Could not create foreach backend object with a specification of class '", class(object)[1L], "'")
	}
)

formatDoName <- function(x){
	
	# numeric values are resolved as doParallel
	if( is.numeric(x) ) x <- 'PAR'
	if( is.character(x) ){
		# use upper case if not already specified as 'do*'
		if( !grepl("^do", x) ){
			x <- toupper(x)
			# special treatment for doParallel
			if( x %in% c('PAR', 'PARALLEL') ) x <- 'Parallel'
		}
		# stick prefix 'do' (removing leading 'do' if necessary)
		str_c('do', sub('^do', '', x))
	}else 
		''
}
#' Creates a foreach backend object based on its name.
setMethod('ForeachBackend', 'character', 
	function(object, ...){
		
		object <- formatDoName(object)
		
		# build S3 class name
		s3class <- str_c(object, "_backend")
		
		# create empty S3 object
		obj <- structure(list(name=object, data=list(...))
						, class=c(s3class, 'foreach_backend'))

		# give a chance to a backend-specific ForeachBackend factory method
		# => this will generally fill the object with the elements suitable
		# to be used in a call to foreach::setDoPar: fun, data, info
		# and possibly change the name or the object class, e.g. to allow 
		# subsequent argument-dependent dispatch.
		obj <- ForeachBackend(obj, ...)
		
		# check the registration routine is available
		.foreach_regfun(obj$name)
		
		# set data slot if not already set by the backend-specific method
		if( is.null(obj$data) || (length(obj$data) == 0L && nargs()>1L) ) 
			obj$data <- list(...)
		
		# return object
		obj
	}
)
#' Creates a foreach backend object for the currently registered backend.
setMethod('ForeachBackend', 'missing', 
	function(object, ...){
		be <- getDoParName()
		data <- getDoBackend()
		bdata <- data$data
		res <- if( !is.null(bdata) ) do.call(ForeachBackend, c(list(be, bdata), ...))
		else ForeachBackend(be, ...)
		if( !is.null(data$cleanup) ) res$cleanup <- data$cleanup
		res
	}
)
#' Dummy method that returns \code{NULL}, defined for correct dispatch.
setMethod('ForeachBackend', 'NULL', function(object, ...){ NULL })

setOldClass('cluster')
#' Creates a doParallel foreach backend that uses the cluster described in 
#' \code{object}.
setMethod('ForeachBackend', 'cluster', 
	function(object, ...){
		ForeachBackend('doParallel', cl=object)
	}
)
#' Creates a doParallel foreach backend with \code{object} processes.
setMethod('ForeachBackend', 'numeric', 
	function(object, ...){
		# check numeric specification
		if( length(object) == 0L )
			stop("invalid number of cores specified as a backend [empty]")
		object <- object[1]
		if( object <= 0 )
			stop("invalid negative number of cores [", object, "] specified for backend 'doParallel'")
		
		ForeachBackend('doParallel', cl=object, ...)
	}
)

setOldClass('doParallel_backend')

setMethod('ForeachBackend', 'doParallel_backend',
	function(object, cl, type=NULL){
				
		# set type of cluster if explicitly provided
		if( !is.null(type) ) object$data$type <- type
        
        # required registration data
		# NB: a function doParallel:::doParallel should exist and do the same 
		# thing as parallel::registerDoParallel without registering the backend
		#object$fun <- doParallel:::doParallel
#		object$info <- doParallel:::info
        # doParallel:::info has been removed from doParallel since version 1.0.7
		# Reported in Issue #7
        object$info <- getDoParallelInfo(object)
		
		# return object
		object
	}
)

setOldClass('doParallelMC_backend')
#' doParallel-specific backend factory for multicore (fork) clusters
#' 
#' This method is needed since version 1.0.7 of \pkg{doParallel}, which removed 
#' internal function \code{info} and defined separate backend names for mc and snow clusters.
setMethod('ForeachBackend', 'doParallelMC_backend',
    function(object, ...){
	
        object$info <- getDoParallelInfo('mc')
        object$name <- 'doParallel'
        # return object
        object
    }
)

setOldClass('doParallelSNOW_backend')
#' doParallel-specific backend factory for SNOW clusters.
#' 
#' This method is needed since version 1.0.7 of \pkg{doParallel}, which removed 
#' internal function \code{info} and defined separate backend names for mc and snow clusters.
setMethod('ForeachBackend', 'doParallelSNOW_backend',
    function(object, ...){
        
        object$info <- getDoParallelInfo('snow')
        object$name <- 'doParallel'
        # return object
        object
    }
)

getDoParallelType <- function(x){
    
    
    cl <- x$data[['cl']]
    if( is.null(cl) && length(x$data) && (is.null(names(x$data)) || names(x$data)[[1L]] == '') )
        cl <- x$data[[1L]]
    if ( is.null(cl) || is.numeric(cl) ) {
        if (.Platform$OS.type == "windows" || (!is.null(x$data$type) && !identical(x$data$type, 'FORK')) ) 'snow'
        else 'mc'
    }
    else 'snow'
    
}

getDoParallelInfo <- function(x, ...){
    t <- if( isString(x) ) x else getDoParallelType(x, ...)
#    str(t)
    ns <- asNamespace('doParallel')
    if( t == 'mc' ) get('mcinfo', ns)
    else get('snowinfo', ns)
}



setOldClass('doPSOCK_backend')
#' doSNOW-specific backend factory
setMethod('ForeachBackend', 'doPSOCK_backend',
		function(object, cl){
			
			# use all available cores if not otherwise specified
			if( missing(cl) ) cl <- getMaxCores()
			
			# return equivalent doParallel object
			ForeachBackend('doParallel', cl, type='PSOCK')
		}
)

.cl_cleanup <- function(gvar, envir=.GlobalEnv){
	if( !exists(gvar, envir = envir) ) return()
	cl <- get(gvar, envir = envir)
	try( parallel::stopCluster(cl), silent=TRUE)
	rm(list=gvar, envir = envir)
	TRUE
} 

cleanupCluster <- function(x, cl, stopFun=NULL){
	
	function(){
		
		if( is(x, 'doParallel_backend') ){
			
            # On non-Windows machines registerDoParallel(numeric) will use 
            # parallel::mclapply with `object` cores (no cleanup required).
            # On Windows doParallel::registerDoParallel(numeric) will create a 
            # SOCKcluster with `object` cores.
            # => Windows needs a cleanup function that will stop the cluster 
            # when another backend is registered.
            # Fortunately doParallel::registerDoParallel assign the cluster object 
            # to the global variable `.revoDoParCluster`
            if( .Platform$OS.type == "windows" ){
				.cl_cleanup(".revoDoParCluster")
			}
		}
		
		if( is.null(stopFun) ) stopFun <- parallel::stopCluster 
		# stop cluster
		stopFun(cl)
		TRUE
	}
}

#' @export
register.doParallel_backend <- function(x, ...){
	
	# start cluster if numeric specification and type is defined
	cl <- x$data[[1]]
  if( is.numeric(cl) && (.Platform$OS.type == 'windows' || !is.null(x$data$type)) ){
		names(x$data)[1L] <- 'spec'
		# start cluster
		clObj <- do.call(parallel::makeCluster, x$data)
		x$data <- list(clObj)
		# setup cleanup procedure
		x$cleanup <- cleanupCluster(x, clObj)
	}
	# register
	register.foreach_backend(x, ...)
}



isMPIBackend <- function(x, ...){
	b <- if( missing(x) ) ForeachBackend(...) else ForeachBackend(object=x, ...)
	if( is.null(b) ) FALSE
	else if( identical(b$name, 'doMPI') ) TRUE 
	else if( length(b$data) ){
		is(b$data[[1]], 'MPIcluster') || is(b$data[[1]], 'mpicluster')
	}else FALSE
}

#' @export
register.doMPI_backend <- function(x, ...){
	
	if( length(x$data) && isNumber(cl <- x$data[[1]]) ){
		clObj <- doMPI::startMPIcluster(cl)
		x$data[[1]] <- clObj
		# setup cleanup procedure
		x$cleanup <- cleanupCluster(x, clObj, doMPI::closeCluster)
	}
	# register
	register.foreach_backend(x, ...)
}

setOldClass('mpicluster')
#' Creates a doMPI foreach backend that uses the MPI cluster described in 
#' \code{object}.
setMethod('ForeachBackend', 'mpicluster', 
	function(object, ...){
		ForeachBackend('doMPI', cl=object)
	}
)

setOldClass('doMPI_backend')
#' doMPI-specific backend factory
setMethod('ForeachBackend', 'doMPI_backend',
	function(object, cl){
		
		# use all available cores if not otherwise specified
		if( missing(cl) ) cl <- getMaxCores()
				
		# required registration data
		object$fun <- doMPI:::doMPI
		object$info <- doMPI:::info
		
		# return object
		object
	}
)


is.backend <- function(x) is(x, 'foreach_backend')

#' @export
print.foreach_backend <- function(x, ...){
	cat("<foreach backend:", x$name, ">\n", sep='')
	if( length(x$data) ){
		cat("Specifications:\n")
		str(x$data)
	}
}

.foreach_regfun <- function(name){
	
	# early exit for doSEQ
	if( name == 'doSEQ' ) return( registerDoSEQ )
    
	# build name of registration function
	s <- str_c(toupper(substring(name, 1,1)), substring(name, 2))
	funname <- str_c('register', s)
	s3class <- str_c(name, "_backend")
	
	# require definition package
	if( !require.quiet(name, character.only=TRUE) )
		stop("could not find package for foreach backend '", name, "'")
	# check for registering function or generic
	if( is.null(regfun <- getFunction(funname, mustFind=FALSE, where=asNamespace(name))) ){
		if( is.null(regfun <- getS3method('register', s3class, optional=TRUE)) )
			stop("could not find registration routine for foreach backend '", name, "'")
		#			stop("backend '", name,"' is not supported: function "
		#							,"`", regfun, "` and S3 method `register.", s3class, "` not found.")
	}
	regfun
}



setGeneric('getDoParHosts', function(object, ...) standardGeneric('getDoParHosts'))
setOldClass('foreach_backend')

setMethod('getDoParHosts', 'ANY',
	function(object, ...){
		
		be <- if( missing(object) ) ForeachBackend(...) else ForeachBackend(object, ...)
		if( existsMethod('getDoParHosts', class(be)[1L]) ) return( callGeneric(object) )
		
		# default behaviour
		nodename <- setNames(Sys.info()['nodename'], NULL)
			
		if( is.null(be) || is.null(be$data) ) return( NULL )
		# doSEQ
		if( be$name == 'doSEQ' ) 
			return( nodename )
		if( isNumber(be$data) ) 
			return( rep(nodename, be$data) )
		if( length(be$data) && isNumber(be$data[[1]]) ) 
			return( rep(nodename, be$data[[1]]) )
		if( length(be$data) && be$name == 'doParallel' ) 
			return( sapply(be$data[[1L]], '[[', 'host') )
		
		if( !missing(object) ){ # backend passed: register temporarly
			ob <- getDoBackend()
			on.exit( setDoBackend(ob) )
			registerDoBackend(be)
		}
		setNames(unlist(times(getDoParWorkers()) %dopar% { Sys.info()['nodename'] }), NULL)
	}
)


getDoParNHosts <- function(object){
	if( missing(object) ) foreach::getDoParWorkers()
	else{
		length(getDoParHosts(object))
	}
}



setupBackend <- function(spec, backend, optional=FALSE, verbose=FALSE){

	pbackend <- backend
	str_backend <- quick_str(pbackend)
	# early exit: FALSE specification or NA backend means not using foreach at all
	if( isFALSE(spec) || is_NA(pbackend) ) return(FALSE)
	# use doParallel with number of cores if specified in backend
	if( is.numeric(pbackend) ){
		spec <- pbackend
		pbackend <- 'PAR'
	}
	# identify doSEQ calls
	doSEQ <- formatDoName(pbackend) == 'doSEQ'
		
	# custom error function
	pcomp <- is.numeric(spec) && !identical(spec[1], 1)
	errorFun <- function(value=FALSE, stop=FALSE, level=1){
		function(e, ...){
			if( !is(e, 'error') ) e <- list(message=str_c(e, ...))
			
			pref <- if( pcomp ) "Parallel" else "Foreach"
			if( !optional || stop ){
				if( verbose >= level ) message('ERROR')
				stop(pref, " computation aborted: ", e$message, call.=FALSE)
			}else if( verbose >= level ){
				message('NOTE')
				message("# NOTE: ", pref, " computation disabled: ", e$message)
			}
			value
		}
	}
	
	# check current backend if backend is NULL
	if( is.null(pbackend) ){
		if( verbose > 1 ){
			message("# Using current backend ... ", appendLF=FALSE)
		}
		ok <- tryCatch({
			if( is.null(parname <- getDoParName()) )
				stop("argument '.pbackend' is NULL but there is no registered backend")
			if( verbose > 1 ) message('OK [', parname, ']')
			TRUE
		}, error = errorFun())
		if( !ok ) return(FALSE)
		# exit now since there is nothing to setup, nothing should change
		# return NULL so that the backend is not restored on.exit of the parent call.
		return(NA)
	}
	##
	
	# test if requested number of cores is actually available
	NCORES <- getMaxCores(limit=FALSE)
	if( verbose > 2 ) message("# Check available cores ... [", NCORES, ']')
	if( verbose > 2 ) message("# Check requested cores ... ", appendLF=FALSE)
	ncores <- if( doSEQ ) 1L
		else{
			ncores <- tryCatch({
					if( is.numeric(spec) ){
						if( length(spec) == 0L )
							stop("no number of cores specified for backend '", str_backend, "'")
						spec <- spec[1]
						if( spec <= 0L )
							stop("invalid negative number of cores [", spec, "] specified for backend '", str_backend, "'")
						spec
					}else # by default use the 'cores' option or half the number of cores
						getMaxCores() #getOption('cores', ceiling(NCORES/2))
				}, error = errorFun(stop=TRUE))
			if( isFALSE(ncores) ) return(FALSE)
			ncores
		}
	if( verbose > 2 ) message('[', ncores, ']')
	
	# create backend object
	if( verbose > 2 ) message("# Loading backend for specification `", str_backend, "` ... ", appendLF=FALSE)
	newBackend <- tryCatch({
			# NB: limit to the number of cores available on the host 
			if( !doSEQ ) ForeachBackend(pbackend, min(ncores, NCORES))
			else ForeachBackend(pbackend)
		}, error = errorFun(level=3))
	if( isFALSE(newBackend) ) return(FALSE)
	if( verbose > 2 ) message('OK')
	
	if( verbose > 1 ) message("# Check host compatibility ... ", appendLF=FALSE)
	ok <- tryCatch({
		# check if we're not running on MAC from GUI
		if( is.Mac(check.gui=TRUE) && (newBackend$name == 'doMC' || (newBackend$name == 'doParallel' && is.numeric(newBackend$data[[1]]))) ){
			# error only if the parallel computation was explicitly asked by the user
			stop("multicore parallel computations are not safe from R.app on Mac OS X."
					, "\n  -> Use a terminal session, starting R from the command line.")
		}
		TRUE
		}, error = errorFun())
	if( !ok ) return(FALSE)
	if( verbose > 1 ) message('OK')
	
	if( verbose > 1 ) message("# Registering backend `", newBackend$name, "` ... ", appendLF=FALSE)
	# try registering the backend
	oldBackend <- getDoBackend()
	# setup retoration of backend in case of an error
	# NB: the new backend cleanup will happens only 
	# if regsitration succeeds, since the cleanup routine is 
	# setup after the registration by the suitable register S3 method. 
	on.exit( setDoBackend(oldBackend, cleanup=TRUE) )
	
	ov <- lverbose(verbose)
	ok <- tryCatch({
			registerDoBackend(newBackend)
			TRUE
		}
		, error ={
			lverbose(ov)
			errorFun()
		})
	lverbose(ov)
	if( !ok ) return(FALSE)
	if( verbose > 1 ) message('OK')
	
	# check allocated cores if not doSEQ backend
	if( newBackend$name != 'doSEQ' ){
		# test allocated number of cores
		if( verbose > 2 ) message("# Check allocated cores ... ", appendLF=FALSE)
		wcores <- getDoParWorkers()
		if( ncores > 0L && wcores < ncores ){
			if( !optional ){
				errorFun(level=3)("only ", wcores, " core(s) available [requested ", ncores ," core(s)]")
			}else if( verbose > 2 ){
				message('NOTE [', wcores, '/', ncores, ']')
				message("# NOTE: using only ", wcores,
						" core(s) [requested ", ncores ," core(s)]")
			}
		}
		else if( verbose > 2 ){
			message('OK [', wcores, '/', ncores
					, if(ncores != NCORES ) str_c(' out of ', NCORES)
					, ']')
		}
	}
	
	# cancel backend restoration
	on.exit()
	# return old backend
	oldBackend
}


# add extra package bigmemory and synchronicity on Unix platforms
if( .Platform$OS.type != 'windows' ){
	setPackageExtra('install.packages', 'bigmemory', pkgs='bigmemory')
	setPackageExtra('install.packages', 'synchronicity', pkgs='synchronicity')
}
# add new option: shared.memory that indicates if one should try using shared memory
# to speed-up parallel computations.
.OPTIONS$newOptions(shared.memory = (.Platform$OS.type != 'windows' && !is.Mac()))



setupSharedMemory <- function(verbose){
	
	if( verbose > 1 ) message("# Check shared memory capability ... ", appendLF=FALSE)
	# early exit if option shared is off
	if( !KINOMO.getOption('shared.memory') ){
		if( verbose > 1 ) message('SKIP [disabled]')
		return(FALSE)
	}
	# early exit if foreach backend is doMPI: it is not working, not sure why
	if( isMPIBackend() ){
		if( verbose > 1 ) message('SKIP [MPI cluster]')
		return(FALSE)
	}
	# not on Windows
	if( .Platform$OS.type == 'windows' ){
		if( verbose > 1 ) message('SKIP [Windows OS]')
		return(FALSE)
	}
	
	if( !require.quiet('bigmemory', character.only=TRUE) ){
		if( verbose > 1 ){
			message('NO', if( verbose > 2 ) ' [Package `bigmemory` required]')
		}
		return(FALSE)
	}
	if( !require.quiet('synchronicity', character.only=TRUE) ){
		if( verbose > 1 ){
			message('NO', if( verbose > 2 ) ' [Package `synchronicity` required]')
		}
		return(FALSE)
	}
	if( verbose > 1 ) message('YES', if( verbose > 2 ) ' [synchronicity]')
	TRUE
}

is.doSEQ <- function(){
	dn <- getDoParName()
	is.null(dn) || dn == 'doSEQ'
}

#' \code{setupTempDirectory} creates a temporary directory to store the best fits computed on each host.
#' It ensures each worker process has access to it.
#' 
#' @rdname setup
setupTempDirectory <- function(verbose){
	
	# - Create a temporary directory to store the best fits computed on each host
	KINOMO_TMPDIR <- tempfile('KINOMO_', getwd())
	if( verbose > 2 ) message("# Setup temporary directory: '", KINOMO_TMPDIR, "' ... ", appendLF=FALSE)
	dir.create(KINOMO_TMPDIR)
	if( !is.dir(KINOMO_TMPDIR) ){
		if( verbose > 2 ) message('ERROR')
		KINOMO_stop('KINOMO', "could not create temporary result directory '", KINOMO_TMPDIR, "'")
	}
	
	on.exit( unlink(KINOMO_TMPDIR, recursive=TRUE) )
	# ensure that all workers can see the temporary directory
	wd <- times(getDoParWorkers()) %dopar% {
		if( !file_test('-d', KINOMO_TMPDIR) )
			dir.create(KINOMO_TMPDIR, recursive=TRUE)
		file_test('-d', KINOMO_TMPDIR)
	}
	# check it worked
	if( any(!wd) ){
		if( verbose > 2 ) message('ERROR')
		KINOMO_stop('KINOMO', "could not create/see temporary result directory '", KINOMO_TMPDIR, "' on worker nodes ", str_out(which(!wd), Inf))	
	}
	if( verbose > 2 ) message('OK')
	on.exit()
	KINOMO_TMPDIR
}

#' Utilities for Parallel Computations
#'
#' 
#' @rdname parallel
#' @name parallel-KINOMO 
NULL


ts_eval <- function(mutex = synchronicity::boost.mutex(), verbose=FALSE){
	
	
	requireNamespace('bigmemory')
	#library(bigmemory)
	requireNamespace('synchronicity')
	#library(synchronicity)
	# describe mutex if necessary
	.MUTEX_DESC <- 
			if( is(mutex, 'boost.mutex') ) synchronicity::describe(mutex)
			else mutex
	
	loadpkg <- TRUE
	function(expr, envir=parent.frame()){
		
		# load packages once
		if( loadpkg ){
			requireNamespace('bigmemory')
			#library(bigmemory)
			requireNamespace('synchronicity')
			#library(synchronicity)
			loadpkg <<- FALSE
		}
		MUTEX <- synchronicity::attach.mutex(.MUTEX_DESC)
		synchronicity::lock(MUTEX)
		if( verbose )
			message('#', Sys.getpid(), " - START mutex: ", .MUTEX_DESC@description$shared.name)
		ERROR <- "### <Error in mutex expression> ###\n"
		on.exit({
			if( verbose ){
				message(ERROR, '#', Sys.getpid(), " - END mutex: ", .MUTEX_DESC@description$shared.name)
			}
			synchronicity::unlock(MUTEX)
		})
		
		eval(expr, envir=envir)
		
		ERROR <- NULL
	}	
}


ts_tempfile <- function(pattern = "file", ..., host=TRUE, pid=TRUE){
	if( host ) pattern <- c(pattern, Sys.info()['nodename'])
	if( pid ) pattern <- c(pattern, Sys.getpid())
	tempfile(paste(pattern, collapse='_'), ...)
}


hostfile <- function(pattern = "file", tmpdir=tempdir(), fileext='', host=TRUE, pid=TRUE){
	if( host ) pattern <- c(pattern, Sys.info()['nodename'])
	if( pid ) pattern <- c(pattern, Sys.getpid())
	file.path(tmpdir, str_c(paste(pattern, collapse='.'), fileext))
}


gVariable <- function(init, shared=FALSE){
	
	if( shared ){ # use bigmemory shared matrices
		if( !is.matrix(init) )
			init <- as.matrix(init)
		requireNamespace('bigmemory')
		#library(bigmemory)
		DATA <- bigmemory::as.big.matrix(init, type='double', shared=TRUE)
		DATA_DESC <- bigmemory::describe(DATA)
	}else{ # use variables assigned to .GlobalEnv
		DATA_DESC <- basename(tempfile('.gVariable_'))
	}
	
	.VALUE <- NULL
	.loadpkg <- TRUE
	function(value){
		
		# load packages once
		if( shared && .loadpkg ){
			requireNamespace('bigmemory')
			#library(bigmemory)
			.loadpkg <<- FALSE	
		}
		
		# if shared: attach bigmemory matrix from its descriptor object
		if( shared ){
			DATA <- bigmemory::attach.big.matrix(DATA_DESC)
		}
		
		if( missing(value) ){# READ ACCESS
			if( !shared ){
				# initialise on first call if necessary
				if( is.null(.VALUE) ) .VALUE <<- init
				# return variable
				.VALUE
			}else 
				DATA[]
			
		}else{# WRITE ACCESS
			if( !shared ) .VALUE <<- value 
			else DATA[] <- value
			
		}
	}
}

#' \code{setupLibPaths} add the path to the KINOMO package to each workers' libPaths. 
#' 
#' @param pkg package name whose path should be exported the workers.
#' 
#' @rdname setup
setupLibPaths <- function(pkg='KINOMO', verbose=FALSE){
	
	# do nothing in sequential mode
	if( is.doSEQ() ) return( character() )
	
	if( verbose ){
		message("# Setting up libpath on workers for package(s) "
			, str_out(pkg, Inf), ' ... ', appendLF=FALSE)
	}
	p <- path.package(pkg)
	if( is.null(p) ) return()
	
	if( !isDevNamespace(pkg) ){ # not a dev package
		plibs <- dirname(p)
		libs <- times(getDoParWorkers()) %dopar% {
			.libPaths(c(.libPaths(), plibs))
		}
		libs <- unique(unlist(libs))
		if( verbose ){
			message("OK\n# libPaths:\n", paste('  ', libs, collapse="\n"))
		}
		libs
		pkg
	}else if( getDoParName() != 'doParallel' || !isNumber(getDoBackend()$data) ){ 
		# devmode: load the package + depends
		if( verbose ){ message("[devtools::load_all] ", appendLF=FALSE) }
		times(getDoParWorkers()) %dopar% {
			capture.output({
				suppressMessages({
					requireNamespace('devtools')
					#library(devtools)
					requireNamespace('bigmemory')
					#library(bigmemory)
					devtools::load_all(p)
				})
			})
		}
		if( verbose ){ message("OK") }
		c('bigmemory', 'rngtools')
	}
	else if( verbose ){
		message("OK")
	}
}



isRNGseed <- function(x){
	is.numeric(x) || 
			( is.list(x) 
				&& is.null(names(x)) 
				&& all(sapply(x, is.numeric)) )
}


setupRNG <- function(seed, n, verbose=FALSE){
	
	if( verbose == 2 ){
		message("# Setting up RNG ... ", appendLF=FALSE)
		on.exit( if( verbose == 2 ) message("OK") )
	}else if( verbose > 2 ) message("# Setting up RNG ... ")
	
	if( verbose > 3 ){
		message("# ** Original RNG settings:")
		showRNG()
	}
	
	# for multiple runs one always uses RNGstreams
	if( n > 1 ){
		
		# seeding with numeric values only
		if( is.list(seed) && isRNGseed(seed) ){
			if( length(seed) != n )
				stop("Invalid list of RNG seeds: must be of length ", n)
			
			if( verbose > 2 ) message("# Using supplied list of RNG seeds")
			return(seed)
			
		}else if( is.numeric(seed) ){
			
			if( verbose > 2 ){
				message("# Generate RNGStream sequence using seed ("
						, RNGstr(seed), ") ... "
						, appendLF=FALSE)
			}
			res <- RNGseq(n, seed)
			if( verbose > 2 ) message("OK")
			return(res)
			
		}else{ # create a sequence of RNGstream using a random seed
			if( verbose > 2 ){
				message("# Generate RNGStream sequence using a random seed ... "
						, appendLF=FALSE)
			}
			res <- RNGseq(n, NULL)
			if( verbose > 2 ) message("OK")
			return(res)
		}
	}else if( is.numeric(seed) ){ 
		# for single runs: 1-length seeds are used to set the current RNG
		# 6-length seeds are used to set RNGstream
		
		if( !is.vector(seed) ){
			message('ERROR')
			stop("KINOMO::KINOMO - Invalid numeric seed: expects a numeric vector.")
		}
		
		# convert to an integer vector
		seed <- as.integer(seed)
		# immediately setup the RNG in the standard way		
		if( length(seed) == 1L ){
			if( verbose > 2 ){
				message("# RNG setup: standard [seeding current RNG]")
				message("# Seeding current RNG with seed (", seed, ") ... "
						, appendLF=FALSE)
			}
			set.seed(seed)
			if( verbose > 2 ) message("OK")				
			return( getRNG() )
		}else if( length(seed) == 6L ){
			if( verbose > 2 ){
				message("# RNG setup: reproducible [using RNGstream]")
				message("# Generate RNGStream sequence using seed ("
						, RNGstr(seed), ") ... "
						, appendLF=FALSE)
			}
			res <- RNGseq(1, seed)
			setRNG(res)
			if( verbose > 2 ) message("OK")
			return( res )
		}else{
			if( verbose > 2 ){
				message("# RNG setup: directly setting RNG")
				message("# Setting RNG with .Random.seed= ("
						, RNGstr(seed), ") ... "
						, appendLF=FALSE)
			}
			setRNG(seed, verbose > 2)
			if( verbose > 2 ) message("OK")
			return( getRNG() )
		}
		stop("KINOMO::KINOMO - Invalid numeric seed: unexpected error.")
	}else{
		if( verbose > 2 ) message("# RNG setup: standard [using current RNG]")
		NULL
	}
} 


