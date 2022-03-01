
#' @include fixed-terms.R
#' @include KINOMOModel.R
NULL




#' 
setClass('KINOMOfit'
	, representation(
			fit = 'KINOMO', # KINOMO model
			residuals = 'numeric', # residuals from the target matrix
			method = 'character', # method used to compute the factorization
			seed = 'character', # seeding method used to compute the factorization
			rng = 'ANY', # numerical random seed
			distance = '.functionSlotNULL', # method used to compute the distance between the target matrix and its KINOMO estimate
			parameters = 'list', # method used to compute the factorization
			runtime = 'proc_time', # running time to perform the KINOMO
			options = 'list', # run options
			extra = 'list' # extra list of results output by the method
			, call = 'call' # store last call to KINOMO()
	)
	
	, prototype = prototype(
			residuals = numeric(),
			method = '',
			seed = '',
			parameters = list(),
			extra = list()
	)
	
	, validity = function(object){
		
		# slot 'objective' must either be a non-empty character string or a function
		obj <- objective(object)
		if( is.character(obj) && obj == '')
			return(paste("Slot 'objective' must either be a non-empty character string or a function definition", sep=''))
		
		
		# everything went fine: return TRUE
		TRUE
	}
	, contains = 'KINOMO'
)


KINOMOfit <- function(fit=KINOMOModel(), ..., rng=NULL){
				
		# use current RNG settings if not otherwise provided
		if( is.null(rng) )
			rng <- getRNG()
		
		new('KINOMOfit', fit=fit, ..., rng=rng)
}


setMethod('fitted', signature(object='KINOMOfit'),
	function(object, ...){
		fitted(fit(object), ...)
	}
)


setMethod('.basis', signature(object='KINOMOfit'),
	function(object, ...){
		.basis(fit(object), ...)
	}
)

setReplaceMethod('.basis', signature(object='KINOMOfit', value='matrix'), 
	function(object, value){ 
		.basis(fit(object)) <- value
		object
	} 
)


setMethod('.coef', signature(object='KINOMOfit'),
	function(object, ...){
		.coef(fit(object), ...)
	}
)

setReplaceMethod('.coef', signature(object='KINOMOfit', value='matrix'), 
	function(object, value){ 
		.coef(fit(object)) <- value
		object
	} 
)


setMethod('ibterms', 'KINOMOfit', 
	function(object){
		ibterms(fit(object))
	}
)

setMethod('icterms', 'KINOMOfit', 
	function(object){
		icterms(fit(object))
	}
)


#' Returns the offset from the fitted model. 
setMethod('offset', signature(object='KINOMOfit'), 
	function(object){
		offset(fit(object))
	}
)


setMethod('niter', signature(object='KINOMOfit'),
	function(object, ...){
		object@extra$iteration
	}
)

setReplaceMethod('niter', signature(object='KINOMOfit', value='numeric'), 
	function(object, value){
		if( (length(value) != 1) || value < 0  ) 
			stop("KINOMO::niter - invalid value for 'niter': single non-negative value is required.", call.=FALSE) 
		object@extra$iteration <- value
		object
	} 
)

#' Show method for objects of class \code{KINOMOfit}
setMethod('show', 'KINOMOfit', 
	function(object)
	{		
		cat("<Object of class: ", class(object), ">\n", sep='')
		cat(" # Model:\n  ")
		s <- capture.output(show(fit(object)))
		cat(s, sep="\n  ")
		cat(" # Details:\n  ")
		.local <- function(){
			if( algorithm(object) != '' ) cat("algorithm: ", algorithm(object), "\n")
			if( seeding(object) != '' ) cat("seed: ",  seeding(object), "\n")
			
			# initial RNG stream			
			cat("RNG: ", RNGstr(object), "\n", sep='')
	
			# distance/objective function
			svalue <- objective(object)
			svalue <- if( is.function(svalue) ) '<function>' else paste("'", svalue,"'", sep='')
			cat("distance metric: ", svalue, "\n")			
			if( length(residuals(object)) !=0  ) cat("residuals: ",  residuals(object), "\n");
			# show the miscellaneous result values
			if( length(object@misc) > 0L )
				cat("miscellaneous:", str_desc(object@misc, exdent=12L), ". (use 'misc(object)')\n")
			# show the parameters specific to the method		
			if( length(object@parameters) > 0 ){
				cat("parameters:", str_desc(object@parameters, exdent=12L), "\n")
#				p <- sapply(object@parameters, function(x){
#					if( is.vector(x) && length(x) == 1L ) x
#					else paste("<", class(x), ">", sep='')
#				})
#				cat(str_wrap(str_out(p, NA, use.names=TRUE, quote=FALSE), exdent=12), "\n")
			}
			# show number of iterations if present
			if( !is.null(i <- niter(object)) ) cat("Iterations:", i, "\n")
			# show elapsed time if present
			if( length(runtime(object)) > 0 ){ cat("Timing:\n"); show(runtime(object));}
		}
		s <- capture.output(.local())
		cat(s, sep="\n  ")
	}
)




#' @export 
setGeneric('fit', function(object, ...) standardGeneric('fit'))
#' Returns the KINOMO model object stored in slot \code{'fit'}. 
setMethod('fit', 'KINOMOfit', function(object) slot(object, 'fit'))


setGeneric('fit<-', function(object, value) standardGeneric('fit<-'))
#' Updates the KINOMO model object stored in slot \code{'fit'} with a new value.
setReplaceMethod('fit', signature(object='KINOMOfit', value='KINOMO'), 
		function(object, value){ 
			slot(object, 'fit') <- value		
			object # TODO: valid object before returning it (+param check=TRUE or FALSE)
		} 
)


setGeneric('minfit', function(object, ...) standardGeneric('minfit') )
#' Returns the object its self, since there it is the result of a single KINOMO run.
setMethod('minfit', 'KINOMOfit', function(object) object)


#' Returns the type of a fitted KINOMO model.
#' It is a shortcut for \code{modelname(fit(object)}. 
setMethod('modelname', signature(object='KINOMOfit'), 
	function(object)
	{
		modelname(fit(object))
	}
)



#' 
setGeneric('residuals', package='stats')

#' 
setMethod('residuals', 'KINOMOfit', 
	function(object, track=FALSE, niter=NULL, ...){ 
		## IMPORTANT: keep this '...' and do not add a 'method' argument as this
		## one is passed by KINOMOfitX::fit (see bug #159) and is not supposed to be 
		## used
		res <- slot(object, 'residuals')
		if( track ) res 
		else if( is.null(niter) ) tail(res, n=1)
		else res[as.character(niter)]
	} 
)


setGeneric('residuals<-', function(object, ..., value) standardGeneric('residuals<-') )
#' @inline
setReplaceMethod('residuals', 'KINOMOfit',
	function(object, ..., niter=NULL, track=FALSE, value){
		if( track ) slot(object, 'residuals') <- value
		else{
			if( !is.null(niter) ) value <- setNames(value, niter)
			slot(object, 'residuals') <- c(slot(object, 'residuals'), value)
		}
		object
	}
)


hasTrack <- function(object, niter=NULL){
	if( is.null(niter) ) length( slot(object, 'residuals') ) > 1
	else !is.na(slot(object, 'residuals')[as.character(niter)])
}


trackError <- function(object, value, niter, force=FALSE){	
	track <- run.options(object, 'error.track')
	track.interval <- run.options(object, 'track.interval')
	
	if( force || (track && niter %% track.interval == 0) ){
		# add the new value to the error track
		last.iter <- names(residuals(object))
		duplicate <- if( !is.null(last.iter) ) niter == last.iter else FALSE
		if( !duplicate ){
			iter <- if( niter >= 0 ) niter
			residuals(object, niter=iter) <- value
		}
	}
	object
}


setMethod('deviance', 'KINOMOfit',
	function(object, y, method, ...){
		
		if( missing(y) ) setNames(residuals(object), NULL)
		else{
			# if missing retrieve the actual distance measure from the KINOMO object
			if( missing(method) ) method = object@distance
			
			# compute the distance between the target and the fitted KINOMO model
			deviance(fit(object), y, method=method, ...)
		}
	}
)

#' Returns the name of the algorithm that fitted the KINOMO model \code{object}.
setMethod('algorithm', 'KINOMOfit', function(object){ object@method } )
#' @inline
setReplaceMethod('algorithm', 'KINOMOfit',
	function(object, value){
		object@method <- value
		object
	}
)

#' Returns the name of the seeding method that generated the starting point
#' for the KINOMO algorithm that fitted the KINOMO model \code{object}.
setMethod('seeding', 'KINOMOfit', function(object){ object@seed } )
#' @inline
setReplaceMethod('seeding', 'KINOMOfit',
	function(object, value){
		object@seed <- value
		object
	}
)


setMethod('objective', signature(object='KINOMOfit'),
	function(object, y){
		
		# when both x and y are missing then returns slot objective
		if( missing(y) ) return(slot(object, 'distance'))
		
		# return the distance computed using the strategy's objective function
		deviance(fit(object), y, method=slot(object, 'distance'))
		
	}
)
#' @inline
setReplaceMethod('objective', signature(object='KINOMOfit', value='ANY'),
	function(object, value){
		slot(object, 'distance') <- value
		validObject(object)
		object
	}
)

#' Returns the CPU time required to compute a single KINOMO fit.
setMethod('runtime', 'KINOMOfit', 
	function(object, ...){ 
		object@runtime
	}
)

#' Identical to \code{runtime}, since their is a single fit. 
setMethod('runtime.all', 'KINOMOfit', getMethod('runtime', 'KINOMOfit'))

###% Access methods to run options.
setGeneric('run.options', function(object, ...) standardGeneric('run.options') )
setMethod('run.options', 'KINOMOfit', 
	function(object, name){
		if( missing(name) ) object@options
		else object@options[[name]]
	}
)
setGeneric('run.options<-', function(object, ..., value) standardGeneric('run.options<-') )
setReplaceMethod('run.options', 'KINOMOfit', 
	function(object, ..., value){
		
		params <- list(...)
		baseError <- 'Setting KINOMO runtime options: ' 
		if ( length(params) == 0 ){
			if( !is.list(value) ) stop(baseError, 'options must be given as a list')
			object@options <- value
			return(object)
		}
		else if ( length(params) > 1 ) stop(baseError, 'options cannot set more than one option at a time')
		name <- params[[1]]
		if( !is.character(name) ) stop(baseError, 'option name must be given as a character string')
		# check if the option exists
		#if( !is.element(name, names(KINOMO.options.runtime())) ) stop(baseError, "option '", name, "' is not defined.")
		
		object@options[[name]] <- value
		object
	}
)
setGeneric('verbose', function(object, ...) standardGeneric('verbose') )
setMethod('verbose', 'KINOMOfit', 
	function(object){
		return(run.options(object, 'verbose') || KINOMO.getOption('debug'))
	}
)

setGeneric('plot', package='graphics' )

setMethod('plot', signature(x='KINOMOfit', y='missing'),
	function(x, y, skip=-1L, ...){
		
		# retrieve the residuals track
		track <- residuals(x, track=TRUE)
		if( length(track) <= 1 ){
			warning(class(x), ' object has no residuals track')
			return(invisible())
		}
		# skip part of the track
		if( skip == -1L && !is.null(names(track)) ) track <- track[names(track)!='0'] # remove initial residual
		else if( skip > 0 ) track <- track[-(1:skip)]
		
		# set default graphical parameters (those can be overriden by the user)
		params <- .set.list.defaults(list(...)
				, xlab='Iterations'
				, ylab=paste('Objective value ('
							, if( is.character(x@distance) ) x@distance else algorithm(x), ')'
							, sep='' )
				, main=paste("KINOMO Residuals\nMethod: ", algorithm(x), " - Rank: ", nbasis(x), sep='')
				, cex.main = 1
				, col='#5555ff', lwd=1.4, type='l', cex=0.5)
		
		do.call('plot', c(list(names(track), track), params))
		points(names(track), track, type='p', cex=0.6, col=params$col)
	}
)


setMethod('summary', signature(object='KINOMOfit'), 
	function(object, ...){
		
		res <- summary(fit(object), ...)
		
		## IMPORTANT: if adding a summary measure also add it in the sorting 
		## schema of method KINOMOfitX::compare to allow ordering on it
		
		# retreive final residuals
		res <- c(res, residuals=as.numeric(residuals(object)))
		# nb of iterations
		res <- c(res, niter=as.integer(niter(object)) )
		# runtime
		t <- runtime(object)
		utime <- as.numeric(t['user.self'] + t['user.child'])
		res <- c(res, cpu=utime, cpu.all=utime, nrun=1)		
		
		# return result
		return(res)
	}
)

#' Compares two KINOMO models when at least one comes from a KINOMOfit object, 
#' i.e. an object returned by a single run of \code{\link{KINOMO}}. 
setMethod('KINOMO.equal', signature(x='KINOMOfit', y='KINOMO'), 
		function(x, y, ...){
			KINOMO.equal(fit(x), y, ...)
		}
)
#' Compares two KINOMO models when at least one comes from a KINOMOfit object, 
#' i.e. an object returned by a single run of \code{\link{KINOMO}}.
setMethod('KINOMO.equal', signature(x='KINOMO', y='KINOMOfit'), 
		function(x, y, ...){
			KINOMO.equal(x, fit(y), ...)
		}
)
#' Compares two fitted KINOMO models, i.e. objects returned by single runs of 
#' \code{\link{KINOMO}}.
setMethod('KINOMO.equal', signature(x='KINOMOfit', y='KINOMOfit'), 
		function(x, y, ...){
			KINOMO.equal(fit(x), fit(y), ...)
		}
)

