#' @include KINOMOStrategy-class.R
#' @include KINOMOfit-class.R
NULL

# Define union class for generalised function slots, e.g., slot 'KINOMOStrategyIterative::Stop' 
setClassUnion('.GfunctionSlotNULL', c('character', 'integer', 'numeric', 'function', 'NULL'))



#'  
setClass('KINOMOStrategyIterative'
	, representation(
                onInit = '.functionSlotNULL',
				Update = '.functionSlot', # update method	
				Stop = '.GfunctionSlotNULL', # method called just after the update
				onReturn = '.functionSlotNULL' # method called just before returning the resulting KINOMO object
				)	
  , prototype=prototype(
          		onInit = NULL
				, Update = ''
				, Stop = NULL
				, onReturn = NULL
			)
	, contains = 'KINOMOStrategy'
	, validity = function(object){
		
		if( is.character(object@Update) && object@Update == '' )
			return("Slot 'Update' is required")
		
		# check the arguments of methods 'Update' and 'Stop'
		# (except for the 3 mandatory ones)
		n.update <- names(formals(object@Update))
		
		# at least 3 arguments for 'Update'
		if( length(n.update) < 3 ){
			return(str_c("Invalid 'Update' method - must have at least 3 arguments: ",
						"current iteration number [i], ",
						"target matrix [y], ",
						"current KINOMO model iterate [x]"))
		}
		
		n.update <- n.update[-seq(3)]
		# argument '...' must be present in method 'Update'
		if( !is.element('...', n.update) )
			return("Invalid 'Update' method: must have argument '...' (even if not used)")

		# at least 3 arguments for 'Stop'
		if( !is.null(object@Stop) ){
			
			# retrieve the stopping criterion and check its intrinsic validity
			.stop <- tryCatch( KINOMOStop(object@Stop, check=TRUE), 
					error = function(e) return(message(e)))
			
			# Update and Stop methods cannot have overlapping arguments
			n.stop <- names(formals(.stop))
			overlap <- intersect(n.update, n.stop)
			overlap <- overlap[which(overlap!='...')]
			if( length(overlap) > 0 ){
				return(str_c("Invalid 'Update' and 'Stop' methods: conflicting arguments ",
						str_out(overlap, Inf)))
			}
		}
		
		TRUE
	}
)


#' Show method for objects of class \code{KINOMOStrategyIterative}
#' @export
setMethod('show', 'KINOMOStrategyIterative',
	function(object){
		
		#cat('<object of class: KINOMOStrategyIterative>')
		callNextMethod()
		cat(" <Iterative schema>\n")
		# go through the slots
		s.list <- names(getSlots('KINOMOStrategyIterative'))
		s.list <- setdiff(s.list, names(getSlots('KINOMOStrategy')))
		#s.list <- s.list[s.list=='ANY']
#		s.list <- c('Update', 'Stop', 'WrapKINOMO')
		out <-
		sapply(s.list, function(sname){
					svalue <- slot(object,sname)
					svalue <- 
					if( is.function(svalue) ) {
						str_args(svalue, exdent=12)
					} else if( is.null(svalue) ){
						'none'
					} else { 
						paste("'", svalue,"'", sep='')
					}
					str_c(sname, ": ", svalue)
				})
		cat(str_c('  ', out, collapse='\n'), "\n", sep='')
		return(invisible())
	}
)
###% This class is an auxiliary class that defines the strategy's methods by directly callable functions. 
setClass('KINOMOStrategyIterativeX'
	, contains = 'KINOMOStrategyIterative'
	, representation = representation(
				workspace = 'environment' # workspace to use persistent variables accross methods
				)
)


###% Creates a KINOMOStrategyIterativeX object from a KINOMOStrategyIterative object.
xifyStrategy <- function(strategy, workspace=new.env(emptyenv())){	
	
	# do nothing if already executable
	if( is(strategy, 'KINOMOStrategyIterativeX') ) return(strategy)
	
	# first check the strategy's validity
	if( is.character(err <- validObject(strategy, test=TRUE)) ){
		stop("Invalid strategy definition:\n\t- ", err)
	}
	
	# intanciate the KINOMOStrategyIterativeX, creating the strategy's workspace
	strategyX <- new('KINOMOStrategyIterativeX', strategy, workspace=workspace)
	
	# define auxiliary function to preload the 'function' slots in class KINOMOStrategyIterativeX
	preload.slot <- function(strategy, sname, default){
		
		# get the content of the slot
		svalue <- slot(strategy,sname)
		
		# if the slot is valid (i.e. it's a non-empty character string), then process the name into a valid function
		fun <-
		if( is.null(svalue) && !missing(default) ) default
		else if( sname == 'Stop' ) KINOMOStop(svalue)
		else if( is.character(svalue) && nchar(svalue) > 0 ){
			# set the slot with the executable version of the function			
			getFunction(svalue)
		}else if( is.function(svalue) )	svalue		
		else
			stop("KINOMOStrategyIterativeX - could not pre-load slot '", sname, "'")		

		# return the loaded function
		fun
	}
	
	# preload the function slots
	slot(strategyX, 'Update') <- preload.slot(strategyX, 'Update')
	slot(strategyX, 'Stop') <- preload.slot(strategyX, 'Stop', function(strategy, i, target, data, ...){FALSE})
	slot(strategyX, 'onReturn') <- preload.slot(strategyX, 'onReturn', identity)
	
	# load the objective function
	objective(strategyX) <- KINOMODistance(objective(strategy))

	# valid the preloaded object
	validObject(strategyX)
	
	# return the executable strategy 
	strategyX
}


staticVar <- local({
	
	.Workspace <- NULL
	function(name, value, init=FALSE){	
		
		# return last workspace
		if( missing(name) ) return(.Workspace)			
		else if( is.null(name) ){ # reset workspace
			.Workspace <<- NULL
			return()
		} else if( is.environment(name) ){ # setup up static environment			
			KINOMO.debug('Strategy Workspace', "initialize static workspace: ", capture.output(.Workspace), "=", capture.output(name))
			.Workspace <<- name
		}else if( isString(name) && is.environment(.Workspace) ){
			if( missing(value) ){
				get(name, envir=.Workspace, inherits=FALSE)
			}else{
				if( !init || !exists(name, envir=.Workspace, inherits=FALSE) )
				{
					if( init ) KINOMO.debug('Strategy Workspace', "initialize variable '", name, "'")
					assign(name, value, envir=.Workspace)
				}
				# return current value
				get(name, envir=.Workspace, inherits=FALSE)
			}
		}else{
			stop("Invalid KINOMO workspace query: .Workspace=", class(.Workspace), '| name=', name
				, if( !missing(value) ) paste0(' | value=', class(value)))
		}
		
	}
})


setMethod('run', signature(object='KINOMOStrategyIterative', y='matrix', x='KINOMOfit'),
	function(object, y, x, .stop=NULL, maxIter = KINOMO.getOption('maxIter') %||% 2000L, ...){
	
	method <- object
	# override the stop method on runtime
	if( !is.null(.stop) ){
		method@Stop <- KINOMOStop(.stop)
		# honour maxIter unless .stop is an integer and maxIter is not passed
		# either directly or from initial call
		# NB: maxIter may be not missing in the call to run() due to the application 
		# of default arguments from the Strategy within KINOMO(), in which case one does not 
		# want to honour it, since it is effectively missing in the original call. 
		if( is.integer(.stop) && (missing(maxIter) || !('maxIter' %in% names(x@call))) )
			maxIter <- .stop[1]
	}
	
	# debug object in debug mode
	if( KINOMO.getOption('debug') ) show(method)		
	
	#Vc# Define local workspace for static variables
	# this function can be called in the methods to get/set/initialize 
	# variables that are persistent within the strategy's workspace
	.Workspace <- new.env()	
	staticVar(.Workspace)
	on.exit( staticVar(NULL) )
		
	# runtime resolution of the strategy's functions by their names if necessary
	strategyX = xifyStrategy(method, .Workspace)
	run(strategyX, y, x, maxIter=maxIter, ...)
})

#' @rdname KINOMOStrategy
setMethod('run', signature(object='KINOMOStrategyIterativeX', y='matrix', x='KINOMOfit'),
	function(object, y, x, maxIter, ...){
				
	strategy <- object
	v <- y
	seed <- x
	#V!# KINOMOStrategyIterativeX::run
	
	#Vc# Define workspace accessor function
	# this function can be called in the methods to get/set/initialize 
	# variables that are persistent within the strategy's workspace

	
	#Vc# initialize the strategy
	# check validity of arguments if possible
	method.args <- KINOMOFormals(strategy, runtime=TRUE)
	internal.args <- method.args$internals
	expected.args <- method.args$defaults
	passed.args <- names(list(...))
	forbidden.args <- is.element(passed.args, c(internal.args))
	if( any(forbidden.args) ){
		stop("KINOMO::run - Update/Stop method : formal argument(s) "
			, paste( paste("'", passed.args[forbidden.args],"'", sep=''), collapse=', ')
			, " already set internally.", call.=FALSE)
	}
	# !is.element('...', expected.args) && 
	if( any(t <- !pmatch(passed.args, names(expected.args), nomatch=FALSE)) ){
		stop("KINOMO::run - onInit/Update/Stop method for algorithm '", name(strategy),"': unused argument(s) "
			, paste( paste("'", passed.args[t],"'", sep=''), collapse=', '), call.=FALSE)
	}
	# check for required arguments
	required.args <- sapply(expected.args, function(x){ x <- as.character(x); length(x) == 1 && nchar(x) == 0 } )
	required.args <- names(expected.args[required.args])
	required.args <- required.args[required.args!='...']
	
	if( any(t <- !pmatch(required.args, passed.args, nomatch=FALSE)) )
		stop("KINOMO::run - Update/Stop method for algorithm '", name(strategy),"': missing required argument(s) "
			, paste( paste("'", required.args[t],"'", sep=''), collapse=', '), call.=FALSE)
	
	#Vc# Start iterations
	KINOMOData <- seed
	# cache verbose level
	verbose <- verbose(KINOMOData)
	
	# clone the object to allow the updates to work in place
	if( verbose > 1L ) 
		message("# Cloning KINOMO model seed ... ", appendLF=FALSE)
	KINOMOFit <- clone(fit(KINOMOData))
	if( verbose > 1L )
		message("[", C.ptr(fit(KINOMOData)), " -> ", C.ptr(KINOMOFit), "]")		
	
	## onInit
	if( is.function(strategy@onInit) ){
		if( verbose > 1L )	message("# Step 1 - onInit ... ", appendLF=TRUE)
		KINOMOFit <- strategy@onInit(strategy, v, KINOMOFit, ...)
		if( verbose > 1L )	message("OK")
	}
	##
	
	# pre-load slots
	updateFun <- strategy@Update
	stopFun <- strategy@Stop
	
	showNIter.step <- 50L
	showNIter <- verbose && maxIter >= showNIter.step
	if( showNIter ){
		ndIter <- nchar(as.character(maxIter))
		itMsg <- paste0('Iterations: %', ndIter, 'i', "/", maxIter)
		cat(itMsgBck <- sprintf(itMsg, 0))
		itMsgBck <- nchar(itMsgBck)
	}
	i <- 0L
	while( TRUE ){
		
		#Vc# Stopping criteria
		# check convergence (generally do not stop for i=0L, but only initialise static variables
		stop.signal <- stopFun(strategy, i, v, KINOMOFit, ...)
		
		# if the strategy ask for stopping, then stop the iteration
		if( stop.signal || i >= maxIter ) break;
		
		# increment i
		i <- i+1L
		
		if( showNIter && (i==1L || i %% showNIter.step == 0L) ){
			cat(paste0(rep("\r", itMsgBck), sprintf(itMsg, i)))
		}
		
		#Vc# update the matrices
		KINOMOFit <- updateFun(i, v, KINOMOFit, ...)
		
		# every now and then track the error if required
		KINOMOData <- trackError(KINOMOData, deviance(strategy, KINOMOFit, v, ...), niter=i)
				
	}
	if( showNIter ){
		ended <- if( stop.signal ) 'converged' else 'stopped' 
		cat("\nDONE (", ended, " at ",i,'/', maxIter," iterations)\n", sep='')
	}
	
	# force to compute last error if not already done
	KINOMOData <- trackError(KINOMOData, deviance(strategy, KINOMOFit, v, ...), niter=i, force=TRUE)
	
	# store the fitted model
	fit(KINOMOData) <- KINOMOFit
	
	#Vc# wrap up
	# let the strategy build the result
	KINOMOData <- strategy@onReturn(KINOMOData)
	if( !inherits(KINOMOData, 'KINOMOfit') ){
		stop('KINOMOStrategyIterative[', name(strategy), ']::onReturn did not return a "KINOMO" instance [returned: "', class(KINOMOData), '"]')
	}
	
	# set the number of iterations performed
	niter(KINOMOData) <- i
	
	#return the result
	KINOMO.debug('KINOMOStrategyIterativeX::run', 'Done')
	invisible(KINOMOData)
})


#' @export
KINOMOFormals.KINOMOStrategyIterative <- function(x, runtime=FALSE, ...){
	
	strategy <- xifyStrategy(x)
	# from run method
	m <- getMethod('run', signature(object='KINOMOStrategyIterative', y='matrix', x='KINOMOfit'))
	run.args <- allFormals(m)[-(1:3)]
	# onInit
	init.args <- if( is.function(strategy@onInit) ) formals(strategy@onInit)
	# Update
	update.args <- formals(strategy@Update)
	# Stop
	stop.args <- formals(strategy@Stop)
	# spplit internals and 
	internal.args <- names(c(init.args[1:3], update.args[1:3], stop.args[1:4]))
	expected.args <- c(init.args[-(1:3)], update.args[-(1:3)], stop.args[-(1:4)])
	
	if( runtime ){
		# prepend registered default arguments
		expected.args <- expand_list(strategy@defaults, expected.args)
		list(internal=internal.args, defaults=expected.args)
	}else{
		args <- c(run.args, expected.args)
		# prepend registered default arguments
		expand_list(strategy@defaults, args)
	}
}



KINOMO_update.KL.h <- std.divergence.update.h <- function(v, w, h, nbterms=0L, ncterms=0L, copy=TRUE)
{	
	.Call("divergence_update_H", v, w, h, nbterms, ncterms, copy, PACKAGE='KINOMO')
}

KINOMO_update.KL.h_R <- R_std.divergence.update.h <- function(v, w, h, wh=NULL)
{	
	# compute WH if necessary	
	if( is.null(wh) ) wh <- w %*% h
	
	# divergence-reducing KINOMO iterations
	# H_au = H_au ( sum_i [ W_ia V_iu / (WH)_iu ] ) / ( sum_k W_ka ) -> each row of H is divided by a the corresponding colSum of W
	h * crossprod(w, v / wh) / colSums(w)	
}


KINOMO_update.KL.w <- std.divergence.update.w <- function(v, w, h, nbterms=0L, ncterms=0L, copy=TRUE)
{	
	.Call("divergence_update_W", v, w, h, nbterms, ncterms, copy, PACKAGE='KINOMO')
}

KINOMO_update.KL.w_R <- R_std.divergence.update.w <- function(v, w, h, wh=NULL)
{			
	# compute WH if necessary	
	if( is.null(wh) ) wh <- w %*% h
	
	# W_ia = W_ia ( sum_u [H_au A_iu / (WH)_iu ] ) / ( sum_v H_av ) -> each column of W is divided by a the corresponding rowSum of H
	#x2 <- matrix(rep(rowSums(h), nrow(w)), ncol=ncol(w), byrow=TRUE); 
	#w * tcrossprod(v / wh, h) / x2;
	sweep(w * tcrossprod(v / wh, h), 2L, rowSums(h), "/", check.margin = FALSE) # optimize version?
	
}


KINOMO_update.euclidean.h <- std.euclidean.update.h <- 
function(v, w, h, eps=10^-9, nbterms=0L, ncterms=0L, copy=TRUE){
	.Call("euclidean_update_H", v, w, h, eps, nbterms, ncterms, copy, PACKAGE='KINOMO')
}

KINOMO_update.euclidean.h_R <- R_std.euclidean.update.h <- function(v, w, h, wh=NULL, eps=10^-9){
	# compute WH if necessary	
	den <- if( is.null(wh) ) crossprod(w) %*% h
			else{ t(w) %*% wh}
	
	# H_au = H_au (W^T V)_au / (W^T W H)_au
	pmax(h * crossprod(w,v),eps) / (den + eps);
}


KINOMO_update.euclidean.w <- std.euclidean.update.w <-
function(v, w, h, eps=10^-9, nbterms=0L, ncterms=0L, weight=NULL, copy=TRUE){
	.Call("euclidean_update_W", v, w, h, eps, weight, nbterms, ncterms, copy, PACKAGE='KINOMO')
}
#' @rdname KINOMO_update_euclidean
#' @export
KINOMO_update.euclidean.w_R <- R_std.euclidean.update.w <- function(v, w, h, wh=NULL, eps=10^-9){
	# compute WH if necessary	
	den <- if( is.null(wh) ) w %*% tcrossprod(h)
			else{ wh %*% t(h)}
	
	# W_ia = W_ia (V H^T)_ia / (W H H^T)_ia and columns are rescaled after each iteration	
	pmax(w * tcrossprod(v, h), eps) / (den + eps);
}



KINOMOStop <- function(s, check=TRUE){
	
	key <- s
	
	fun <- 
	if( is.integer(key) )	KINOMO.stop.iteration(key)
	else if( is.numeric(key) ) KINOMO.stop.threshold(key)
	else if( is.function(key) ) key
	else if( is.character(key) ){
		# update .stop for back compatibility:
		if( key == 'KINOMO.stop.consensus') key <- 'connectivity'
		
		# first lookup for a `KINOMO.stop.*` function
		key2 <- paste('KINOMO.stop.', key, sep='')
		e <- pkgmaker::packageEnv()
		sfun <- getFunction(key2, mustFind=FALSE, where = e)
		if( is.null(sfun) ) # lookup for the function as such
			sfun <- getFunction(key, mustFind = FALSE, where = e)			
		if( is.null(sfun) )
			stop("Invalid key ['", key,"']: could not find functions '",key2, "' or '", key, "'")
		sfun
	}else if( identical(key, FALSE) ) # create a function that does not stop 
		function(strategy, i, target, data, ...){FALSE}
	else
		stop("Invalid key: should be a function, a character string or a single integer/numeric value. See ?KINOMOStop.")

	# execute if generator (i.e. no arguments)
	if( length(formals(fun)) == 0L ) fun <- fun() 

	# check validity if requested
	if( check ){
		n.stop <- names(formals(fun))
		if( length(n.stop) < 4 ){
			stop("Invalid 'Stop' method - must have at least 4 arguments: ",
					"KINOMO strategy object [strategy], ",
					"current iteration number [i], ",
					"target matrix [y], ",
					"current KINOMO model iterate [x]")
		}
		
		n.stop <- n.stop[-seq(4)]
		# argument '...' must be present in method 'Stop'
		if( !is.element('...', n.stop) )
			stop("Invalid 'Stop' method: must have argument '...' (even if not used)")
	}
	
	# return function
	fun
}


KINOMO.stop.iteration <- function(n){
	
	KINOMO.debug("Using stopping criterion - Fixed number of iterations: ", n)
	if( !is.numeric(n) )
		stop("Invalid argument `n`: must be an integer value")
	if( length(n) > 1 )
		warning("KINOMO::KINOMO - Argument `n` [", deparse(substitute(n)), "] has length > 1: only using the first element.")
	
	.max <- n[1]
	function(object, i, y, x, ...) i >= .max
}


KINOMO.stop.threshold <- function(threshold){	
	
	KINOMO.debug("Using stopping criterion - Stationarity threshold: ", threshold)
	if( !is.numeric(threshold) )
		stop("Invalid argument `threshold`: must be a numeric value")
	if( length(threshold) > 1 )
		warning("KINOMO::KINOMO - Argument `threshold` [", deparse(substitute(threshold)), "] has length > 1: only using the first element.")
	
	eval(parse(text=paste("function(strategy, i, target, data, stationary.th=", threshold, ", ...){
		KINOMO.stop.stationary(strategy, i, target, data, stationary.th=stationary.th, ...)
	}")))
}




KINOMO.stop.stationary <- local({
	
	# static variable
	.last.objective.value <- c(-Inf, Inf)
	.niter <- 0L
	
	.store_value <- function(value){
		.niter <<- .niter + 1L
		.last.objective.value <<- c(max(.last.objective.value[1L], value)
									, min(.last.objective.value[2L], value))
	}
	
	.reset_value <- function(){
		.last.objective.value <<- c(-Inf, Inf)
		.niter <<- 0L
	}
	
	function(object, i, y, x, stationary.th=.Machine$double.eps, check.interval=5*check.niter, check.niter=10L, ...){
		
		# check validity
		if( check.interval < check.niter ){
			stop("Invalid argument values: `check.interval` must always be greater than `check.niter`")
		}
		# initialisation call: compute initial objective value
		if( i == 0L || (i == 1L && is.null(.last.objective.value)) ){
			.reset_value()
			
			# give the chance to update once and estimate from a partial model
			if( is.partial.KINOMO(x) ) return( FALSE )
			
			# compute initial deviance
			current.value <- deviance(object, x, y, ...)
			# check for NaN, i.e. probably infinitely small value (cf. bug reported by Nadine POUKEN SIEWE)
			if( is.nan(current.value) ) return(TRUE)
			
			# store value in static variable for next calls
			.store_value(current.value)
			
			return(FALSE)
		}
		
		# test convergence only every 10 iterations
		if( .niter==0L && i %% check.interval != 0 ) return( FALSE );
		
		# get last objective value from static variable		
		current.value <- deviance(object, x, y, ...)
		# check for NaN, i.e. probably infinitely small value (cf. bug reported by Nadine POUKEN SIEWE)
		if( is.nan(current.value) ) return(TRUE)
		
		# update static variables
		.store_value(current.value)
		
		# once values have been computed for check.niter iterations:
		# check if the difference in the extreme objective values is small enough
		if( .niter == check.niter+1 ){
			crit <- abs(.last.objective.value[1L] - .last.objective.value[2L]) / check.niter
			if( crit <= stationary.th ){
				if( KINOMO.getOption('verbose') ){
					message(crit)
				}
				return( TRUE )
			}
			.reset_value()
		}
		
		# do NOT stop
		FALSE
	}
})


KINOMO.stop.connectivity <- local({
			
	# static variables
	.consold <- NULL
	.inc <- NULL
	
	function(object, i, y, x, stopconv=40, check.interval=10, ...){

		if( i == 0L ){ # initialisation call
			# Initialize consensus variables 
			# => they are static variables within the strategy's workspace so that
			# they are persistent and available throughout across the calls
			p <- ncol(x)
			.consold <<- matrix(0, p, p)
			.inc <<- 0
			return(FALSE)
		}
	
		# test convergence only every 10 iterations
		if( i %% check.interval != 0 ) return( FALSE );
		
		# retrieve metaprofiles
		h <- coef(x, all=FALSE)
			
		# construct connectivity matrix
		index <- apply(h, 2, function(x) which.max(x) )
		cons <- outer(index, index, function(x,y) ifelse(x==y, 1,0));

		changes <- cons != .consold
		if( !any(changes) ) .inc <<- .inc + 1 # connectivity matrix has not changed: increment the count
		else{
			.consold <<- cons;
			.inc <<- 0;                         # else restart counting
		}
		
		# prints number of changing elements 		
		#if( verbose(x) ) cat( sprintf('%d ', sum(changes)) ) 
		#cat( sprintf('%d ', sum(changes)) )
										
		# assume convergence is connectivity stops changing 
		if( .inc > stopconv ) return( TRUE );
		
		# do NOT stop
		FALSE
	}
})

