#' @include KINOMOstd-class.R
#' @include KINOMOSet-class.R
#' @include registry-seed.R
#' @include registry-algorithms.R
#' @include parallel.R

NULL

setGeneric('KINOMO', function(x, rank, method, ...) standardGeneric('KINOMO') )

setMethod('KINOMO', signature(x='data.frame', rank='ANY', method='ANY'), 
	function(x, rank, method, ...)
	{
		# replace missing values by NULL values for correct dispatch
		if( missing(method) ) method <- NULL
		if( missing(rank) ) rank <- NULL
		
		# apply KINOMO to the the data.frame converted into a matrix	
		KINOMO(as.matrix(x), rank, method, ...)
	}
)

setMethod('KINOMO', signature(x='matrix', rank='numeric', method='NULL'), 
		function(x, rank, method, seed=NULL, model=NULL, ...)
		{
			
			# a priori the default method will be used
			method <- KINOMO.getOption('default.algorithm')
			
			# use default seeding method if seed is missing
			if( is.null(seed) ){
#				seed <- KINOMO.getOption('default.seed')
			}else{
				# get reference object from which to infer model type
				refobj <- if( is.KINOMO(seed) ) seed else if( is.KINOMO(model) ) model
				
				if( !is.null(refobj) ){
					mtype <- modelname(refobj)
					# try to find the algorithm suitable for the seed's KINOMO model
					method.potential <- selectKINOMOMethod(model=mtype, exact=TRUE, quiet=TRUE)
					if( is.null(method.potential) )
						stop("KINOMO::KINOMO - Found no algorithm defined for model '", mtype, "'")
					
					if( length(method.potential) == 1 ) # only one to choose
						method <- method.potential
					else if( !is.element(method, method.potential) ){# several options, none is default
						method <- method.potential[1]
						warning("KINOMO::KINOMO - Selected algorithm '", method, "' to fit model '", mtype, "'."
							, "\n  Alternatives are: "
							, str_out(method.potential[-1], Inf)
							, call.=FALSE, immediate.=TRUE)
					}
				}
			}
			
			KINOMO(x, rank, method, seed=seed, model=model, ...)
		}
)



setMethod('KINOMO', signature(x='matrix', rank='numeric', method='list'), 
	function(x, rank, method, ..., .parameters = list())
	{
		# apply each KINOMO algorithm
		k <- 0
		n <- length(method)
        
        # setup/check method specific parameters
        ARGS <- NULL
        .used.parameters <- character()
        if( !is.list(.parameters) )
            stop("KINOMO::KINOMO - Invalid value for argument `.parameters`: must be a named list.")
        if( length(.parameters) && (is.null(names(.parameters)) || any(names(.parameters) == '')) )
            stop("KINOMO::KINOMO - Invalid value for argument `.parameters`: all elements must be named.") 
        
        t <- system.time({
			res <- lapply(method, 
				function(meth, ...){
					k <<- k+1
					methname <- if( isString(meth) ) meth else name(meth)
					cat("Compute KINOMO method '", methname, "' [", k, "/", n, "] ... ", sep='')
					# restore RNG on exit (except after last method)
					# => this ensures the methods use the same stochastic environment
					orng <- RNGseed()
					if( k < n ) on.exit( RNGseed(orng), add = TRUE)
					
                    # look for method-specific arguments
                    i.param <- 0L
                    if( length(.parameters) ){
                        i.param <- charmatch(names(.parameters), methname)
                        if( !length(i.param <- seq_along(.parameters)[!is.na(i.param)]) )
                            i.param <- 0L
                        else if( length(i.param) > 1L ){
                            stop("Method name '", methname, "' matches multiple method-specific parameters "
                                    , "[", str_out(names(.parameters)[i.param], Inf), "]")
                        }
                    }
					#o <- capture.output( 
                        if( !i.param ){
                            res <- try( KINOMO(x, rank, meth, ...) , silent=TRUE)
                        }else{
                            if( is.null(ARGS) ) ARGS <<- list(x, rank, ...)
                            .used.parameters <<- c(.used.parameters, names(.parameters)[i.param])
                            res <- try( do.call(KINOMO, c(ARGS, method = meth, .parameters[[i.param]])) 
                                        , silent=TRUE)
                        } 
					#)
					if( is(res, 'try-error') )
						cat("ERROR\n")
					else 
						cat("OK\n")
					return(res)
				}
				, ...)
		})
		
		# filter out bad results
		ok <- sapply(res, function(x){
					if( is(x, 'KINOMO.rank') ) all(sapply(x$fit, isKINOMOfit))
					else isKINOMOfit(x)
			})
		if( any(!ok) ){ # throw warning if some methods raised an error
			err <- lapply(which(!ok), function(i){ paste("'", method[[i]],"': ", res[[i]], sep='')})
			warning("KINOMO::KINOMO - Incomplete results due to ", sum(!ok), " errors: \n- ", paste(err, collapse="- "), call.=FALSE)
		}
		res <- res[ok]
		# TODO error if ok is empty

        # not-used parameters
        if( length(.used.parameters) != length(.parameters) ){
            warning("KINOMO::KINOMO - Did not use methods-specific parameters ", str_out(setdiff(names(.parameters), .used.parameters), Inf))
        }

		# add names to the result list
		names(res) <- sapply(res, function(x){
					if( is(x, 'KINOMO.rank') ) x <- x$fit[[1]]
					algorithm(x)
				})
				
		# return list as is if surveying multiple ranks 
		if( length(rank) > 1 ) return(res)
		
		# wrap the result in a KINOMOList object
		# DO NOT WRAP anymore here: KINOMOfitX objects are used only for results of multiple runs (single method)
		# the user can still join on the result if he wants to
		#res <- join(res, runtime=t)
		res <- new('KINOMOList', res, runtime=t)
		
		# return result
		return(res)
	}
)



setMethod('KINOMO', signature(x='matrix', rank='numeric', method='character'),
function(x, rank, method, ...)
{	
	# if there is more than one methods then treat the vector as a list
	if( length(method) > 1 ){
		return( KINOMO(x, rank, as.list(method), ...) )
	}
	
	# create the KINOMOStrategy from its name
	strategy <- KINOMOAlgorithm(method)		
	# apply KINOMO using the retrieved strategy		
	KINOMO(x, rank, method=strategy, ...)
}
)


setMethod('KINOMO', signature(x='matrix', rank='numeric', method='function'),
	function(x, rank, method, seed, model='KINOMOstd', ..., name, objective='euclidean', mixed=FALSE){

		model_was_a_list <- is.list(model)
		if( is.character(model) )
			model <- list(model=model) 
		if( !is.list(model) ){
			stop("KINOMO - Invalid argument `model`: must be NULL or a named list of initial values for slots in an KINOMO model.")
		}
				
		
		# arguments passed to the call to KINOMOStrategyFunction
		strat <- list('KINOMOStrategyFunction'
					, algorithm = method
					, objective = objective
					, mixed = mixed[1]
					)
		
		## Determine type of KINOMO model associated with the KINOMOStrategy
		# All elements of `model` (except the model class) will be passed to 
		# argument `model` of the workhorse `KINOMO` method, which will use them  
		# to create the KINOMO model in a call to `KINOMOModel`
		if( length(model) > 0L ){
			if( !is.null(model$model) ){
				strat$model <- model$model
				model$model <- NULL
			}else if( isKINOMOclass(model[[1]]) ){
				strat$model <- model[[1]]
				# use the remaining elements to instanciate the KINOMO model
				model <- model[-1]
			}
			# all elements must be named
			if( !hasNames(model, all=TRUE) ){
				stop("KINOMO::KINOMO - Invalid argument `model`: all elements must be named, except the first one which must then be an KINOMO model class name")
			}
		}
		##
		
		# if name is missing: generate a temporary unique name
		if( missing(name) ) name <- basename(tempfile("KINOMO_"))
		# check that the name is not a registered name
		if( existsKINOMOMethod(name) )
			stop("Invalid name for custom KINOMO algorithm: '",name,"' is already a registered KINOMO algorithm")
		strat$name <- name

		# create KINOMOStrategy
		strategy <- do.call('new', strat)
		# full validation of the strategy
		validObject(strategy, complete=TRUE)
		
		if( missing(seed) ) seed <- NULL
		if( !model_was_a_list && length(model) == 0L ) model <- NULL
		# call method 'KINOMO' with the new object
		KINOMO(x, rank, strategy, seed=seed, model=model, ...)
	}
)


setMethod('KINOMO', signature(x='matrix', rank='KINOMO', method='ANY'),
	function(x, rank, method, seed, ...){
		
		if( !missing(seed) ){		
			if( isNumber(seed) ){
				set.seed(seed)
			}else if( !is.null(seed) ){
				warning("KINOMO::KINOMO - Discarding value of argument `seed`: directly using KINOMO model supplied in `rank` instead.\n"
						, "  If seeding is necessary, please use argument `model` pass initial model slots, which will be filled by the seeding method.")
			}
#			# pass the model via a one-off global variable
#			.KINOMO_InitModel(rank)
		}
		
		# replace missing method by NULL for correct dispatch
		if( missing(method) ) method <- NULL
		
		KINOMO(x, nbasis(rank), method, seed=rank, ...)
	}
)
.KINOMO_InitModel <- oneoffVariable()


setMethod('KINOMO', signature(x='matrix', rank='NULL', method='ANY'),
	function(x, rank, method, seed, ...){
		
		if( missing(seed) || !is.KINOMO(seed) )
			stop("KINOMO::KINOMO - Argument `seed` must be an KINOMO model when argument `rank` is missing.")
		
		# replace missing method by NULL for correct dispatch
		if( missing(method) ) method <- NULL
		
		KINOMO(x, nbasis(seed), method, seed=seed, ...)
	}
)

setMethod('KINOMO', signature(x='matrix', rank='missing', method='ANY'),
	function(x, rank, method, ...){
		# replace missing method by NULL for correct dispatch
		if( missing(method) ) method <- NULL
		KINOMO(x, NULL, method, ...)
	}
)

setMethod('KINOMO', signature(x='matrix', rank='numeric', method='missing'),
	function(x, rank, method, ...){
		KINOMO(x, rank, NULL, ...)
	}
)

setMethod('KINOMO', signature(x='matrix', rank='matrix', method='ANY'), 
	function(x, rank, method, seed, model=list(), ...)
	{
		if( is.character(model) )
			model <- list(model=model) 
		if( !is.list(model) )
			stop("KINOMO - Invalid argument `model`: must be NULL or a named list of initial values for slots in an KINOMO object.")
		if( !hasNames(model, all=TRUE) )
			stop("KINOMO - Invalid argument `model`: all elements must be named")
			
		# remove rank specification if necessary
		if( !is.null(model$rank) ){
			warning("KINOMO - Discarding rank specification in argument `model`: use value inferred from matrix supplied in argument `rank`")
			model$rank <- NULL
		}
		# check compatibility of dimensions
		newseed <- 
		if( nrow(rank) == nrow(x) ){
			# rank is the initial value for the basis vectors
			if( length(model)==0L ) KINOMOModel(W=rank)
			else{
				model$W <- rank
				do.call('KINOMOModel', model)
			}
		}else if( ncol(rank) == ncol(x) ){
            # rank is the initial value for the mixture coefficients
            if( length(model)==0L ) KINOMOModel(H=rank)
            else{
                model$H <- rank
                do.call('KINOMOModel', model)
            }
        }else
			stop("KINOMO - Invalid argument `rank`: matrix dimensions [",str_out(dim(x),sep=' x '),"]"
				, " are incompatible with the target matrix [", str_out(dim(x),sep=' x '),"].\n"
				, "  When `rank` is a matrix it must have the same number of rows or columns as the target matrix `x`.")
		
		# replace missing values by NULL values for correct dispatch
		if( missing(method) ) method <- NULL
		if( missing(seed) ) seed <- NULL
		#KINOMO(x, nbasis(newseed), method, seed=seed, model=newseed, ...)
		KINOMO(x, newseed, method, seed=seed, ...)
	}
)

setMethod('KINOMO', signature(x='matrix', rank='data.frame', method='ANY'),
	function(x, rank, method, ...){
		# replace missing values by NULL values for correct dispatch
		if( missing(method) ) method <- NULL
		
		KINOMO(x, as.matrix(rank), method, ...)
	}
)


setMethod('KINOMO', signature(x='formula', rank='ANY', method='ANY'),
	function(x, rank, method, ..., model=NULL){
		# replace missing values by NULL values for correct dispatch
		if( missing(method) ) method <- NULL
		if( missing(rank) ) rank <- NULL
		
		# if multiple numeric rank: use KINOMORestimateRank
		if( is.vector(rank) && is.numeric(rank)  ){
			if( length(rank) > 1L ){
				return( KINOMOEstimateRank(x, rank, method, ..., model=model) )
			}
		}
		
		# build formula based model
		model <- KINOMOModel(x, rank, data=model)
		KINOMO(attr(model, 'target'), nbasis(model), method, ..., model=model)
	}
)

.as.numeric <- function(x){
	suppressWarnings( as.numeric(x) )
}

.translate.string <- function(string, dict){
	
	res <- list()
	dict <- as.list(dict)
	if( nchar(string) == 0 ) return(res)
	opt.val <- TRUE
	last.key <- NULL
	buffer <- ''
	lapply(strsplit(string, '')[[1]], 
		function(c){
			if( c=='-' ) opt.val <<- FALSE
			else if( c=='+' ) opt.val <<- TRUE
			else if( opt.val && !is.na(.as.numeric(c)) )
				buffer <<- paste(buffer, c, sep='')
			else if( !is.null(dict[[c]]) ){
				# flush the buffer into the last key if necessary
				if( nchar(buffer) > 0 && !is.null(last.key) && !is.na(buffer <- .as.numeric(buffer)) ){
					res[[dict[[last.key]]]] <<- buffer
					buffer <<- ''
				}
					
				res[[dict[[c]]]] <<- opt.val
				last.key <<- c
			}
		}
	)

	# flush the buffer into the last key
	if( nchar(buffer) > 0 && !is.null(last.key) && !is.na(buffer <- .as.numeric(buffer)) )
		res[[dict[[last.key]]]] <- buffer

	# return result	
	return(res)
}


checkErrors <- function(object, element=NULL){
	
	# extract error messages
	errors <- 
			if( is.null(element) ){
				lapply(seq_along(object), function(i){ 
							x <- object[[i]]
							if( is(x, 'error') ) c(i, x)
							else NA 
						})
			}else{
				lapply(seq_along(object), function(i){
							x <- object[[i]][[element, exact=TRUE]]
							if( is(x, 'error') ) c(i, x) 
							else NA
						})
			}
	errors <- errors[!is.na(errors)]
	nerrors <- length(errors)
	res <- list(n = nerrors)
	
	# format messages
	if( nerrors ){
		ierrors <- sapply(errors, '[[', 1L)
		msg <- sapply(errors, '[[', 2L)
		ierrors_unique <- ierrors[!duplicated(msg)]
		res$msg <- str_c("  - ", str_c("run #", ierrors_unique, ': ', msg[ierrors_unique], collapse="\n  - "))
	}
	
	# return error data
	res
}









setMethod('KINOMO', signature(x='matrix', rank='numeric', method='KINOMOStrategy'),
#function(x, rank, method, seed='random', nrun=1, keep.all=FALSE, optimized=TRUE, init='KINOMO', track, verbose, ...)
function(x, rank, method
		, seed=KINOMO.getOption('default.seed'), rng = NULL
		, nrun=if( length(rank) > 1L ) 30 else 1, model=NULL, .options=list()
		, .pbackend=KINOMO.getOption('pbackend')
		, .callback=NULL #callback function called after a run  
		, ...)
{
	fwarning <- function(...) KINOMO_warning('KINOMO', ...)
	fstop <- function(...) KINOMO_stop('KINOMO', ...)
	n <- NULL
	RNGobj <- NULL
	
	# if options are given as a character string, translate it into a list of booleans
	if( is.character(.options) ){
		.options <- .translate.string(.options, 
				c(t='track', v='verbose', d='debug'
				, p='parallel', P='parallel.required'
				, k='keep.all', r='restore.seed', f='dry.run'
				, g='garbage.collect'
				, c='cleanup', S='simplifyCB'
				, R='RNGstream', m='shared.memory'))
	}
	
	# get seeding method from the strategy's defaults if needed
	seed <- defaultArgument(seed, method, KINOMO.getOption('default.seed'), force=is.null(seed))
	.method_defaults <- method@defaults
	.method_defaults$seed <- NULL
	#
	# RNG specification
	if( isRNGseed(seed) ){
		if( !is.null(rng) )
			warning("Discarding RNG specification in argument `rng`: using those passed in argument `seed`.")
		rng <- seed
		seed <- 'random'
	}
	#

	# setup verbosity options
	debug <- if( !is.null(.options$debug) ) .options$debug else KINOMO.getOption('debug')
	verbose <- if( debug ) Inf
				else if( !is.null(.options$verbose) ) .options$verbose
				else KINOMO.getOption('verbose')
	
	# show call in debug mode
	if( debug ){
		.ca <- match.call()
		message('# KINOMO call: ', paste(capture.output(print(.ca)), collapse="\n  "))
	}
	# KINOMO over a range of values: pass the call to KINOMOEstimateRank
	if( length(rank) > 1 ){
		if( verbose <= 1 )
			.options$verbose <- FALSE
		return( KINOMOEstimateRank(x, range = rank, method = method, nrun = nrun
								, seed = seed, rng = rng, model = model
								, .pbackend = .pbackend, .callback = .callback
								, verbose=verbose, .options=.options, ...) )
	}
	
	.OPTIONS <- list()
	# cleanup on exit
	.CLEANUP <- .options$cleanup %||% TRUE
	
	# tracking of objective value
	.OPTIONS$track <- if( !is.null(.options$track) ) .options$track 
					else KINOMO.getOption('track')
	# dry run
	dry.run <- .options$dry.run %||% FALSE 
	# call the garbage collector regularly
	opt.gc <- if( !is.null(.options$garbage.collect) ) .options$garbage.collect
			  else KINOMO.getOption('gc')
	if( is.logical(opt.gc) && opt.gc )
		opt.gc <- ceiling(max(nrun,50) / 3)
	.options$garbage.collect <- opt.gc
	
	# keep results from all runs?
	keep.all <- .options$keep.all %||% FALSE
    # shared memory?
    shared.memory <- if( !is.null(.options$shared.memory) ) .options$shared.memory else KINOMO.getOption('shared.memory')
	# use RNG stream
	.options$RNGstream <- .options$RNGstream %||% TRUE
	
	# discard .callback when not used
	if( is.function(.callback) ){
		w <- if( nrun==1 ) "discarding argument `.callback`: not used when `nrun=1`."
			else if( keep.all )	
				"discarding argument `.callback`: not used when option `keep.all=TRUE`."
		if( !is.null(w) ){
			.callback <- NULL
			fwarning(w, immediate.=TRUE)
		}
		
		# wrap into another function if necessary
		if( is.function(.callback) ){
			# default is to simplify
			.options$simplifyCB <- .options$simplifyCB %||% TRUE 
			args <- formals(.callback)
			if( length(args) <= 2L ){
				if( length(args) < 2L || '...' %in% names(args) ){
					.CALLBACK <- .callback
					.callback <- function(object, i) .CALLBACK(object)
				}
			}
			
			# define post-processing function
			processCallback <- function(res){
				# check errors
				errors <- checkErrors(res, '.callback')
				if( errors$n > 0 ){
					fwarning("All KINOMO fits were successful but ", errors$n, "/", nrun, " callback call(s) threw an error.\n"
							,"# ", if(errors$n>10) "First 10 c" else "C", "allback error(s) thrown:\n"
							, errors$msg
					)
				}
				# add callback values to result list
				sapply(res, '[[', '.callback'
					, simplify=.options$simplifyCB && errors$n == 0L)
			}
		}
	}
	
	## ROLLBACK PROCEDURE
	exitSuccess <- exitCheck()
	on.exit({ 
		if( verbose > 1 ) message("# KINOMO computation exit status ... ", if( exitSuccess() ) 'OK' else 'ERROR')
		if( verbose > 2 ){
			if( exitSuccess() ){
				message('\n## Running normal exit clean up ... ')
			}else{ 
				message('\n## Running rollback clean up ... ')
			}
		}
	}, add=TRUE)
	# RNG restoration on error
	.RNG_ORIGIN <- getRNG()
	on.exit({
		if( !exitSuccess() ){
			if( verbose > 2 ) message("# Restoring RNG settings ... ", appendLF=verbose>3)
			setRNG(.RNG_ORIGIN)
			if( verbose > 3 ) showRNG(indent=' #')
			if( verbose > 2 ) message("OK")
		}
	}, add=TRUE)

	# Set debug/verbosity option just for the time of the run
	old.opt <- KINOMO.options(debug=debug, verbose=verbose, shared.memory = shared.memory);
	on.exit({
		if( verbose > 2 ) message("# Restoring KINOMO options ... ", appendLF=FALSE)
		KINOMO.options(old.opt)
		if( verbose > 2 ) message("OK")
	}, add=TRUE)
	
	# make sure rank is an integer
	rank <- as.integer(rank)
	if( length(rank) != 1 ) fstop("invalid argument 'rank': must be a single numeric value")
	if( rank < 1 ) fstop("invalid argument 'rank': must be greater than 0")
	
	# option 'restore.seed' is deprecated
	if( !is.null(.options$restore.seed) )
		fwarning("Option 'restore.seed' is deprecated and discarded since version 0.5.99.")
	
	if( verbose ){
		if( dry.run ) message("*** fake/dry-run ***")
		message("KINOMO algorithm: '", name(method), "'")
	}
	
	##START_MULTI_RUN
	# if the number of run is more than 1, then call itself recursively
	if( nrun > 1 )
	{
		if( verbose ) message("Multiple runs: ", nrun)
		
		if( verbose > 3 ){
			cat("## OPTIONS:\n")
			sapply(seq_along(.options)
					, function(i){
						r <- i %% 4
						cat(if(r!=1) '\t| ' else "# ", names(.options)[i],': ', .options[[i]], sep='')
						if(r==0) cat("\n# ")
					})
			if( length(.options) %% 4 != 0 )cat("\n")
		}
				
		## OPTIONS: parallel computations 
		# option require-parallel: parallel computation is required if TRUE or numeric != 0
		opt.parallel.required <- !is.null(.options$parallel.required) && .options$parallel.required
		# determine specification for parallel computations 
		opt.parallel.spec <- 
				if( opt.parallel.required ){ # priority over try-parallel 
					# option require-parallel implies and takes precedence over option try-parallel
					.options$parallel.required
				}else if( !is.null(.options$parallel) ) .options$parallel # priority over .pbackend
				else !is_NA(.pbackend) # required only if backend is not trivial

		# determine if one should run in parallel at all: TRUE or numeric != 0, .pbackend not NA
		opt.parallel <- !is_NA(.pbackend) && (isTRUE(opt.parallel.spec) || opt.parallel.spec)
		##
		if( opt.parallel ){
			if( verbose > 1 )
				message("# Setting up requested `foreach` environment: "
						, if( opt.parallel.required ) 'require-parallel' else 'try-parallel'
						, ' [', quick_str(.pbackend) , ']')
			
			
			# switch doMC backend to doParallel
			if( isString(.pbackend, 'MC', ignore.case=TRUE) ){ 
				.pbackend <- 'par'
			}
			# try setting up parallel foreach backend
			oldBackend <- setupBackend(opt.parallel.spec, .pbackend, !opt.parallel.required, verbose=verbose)
			opt.parallel <- !isFALSE(oldBackend)
			# setup backend restoration if using one different from the current one
			if( opt.parallel && !is_NA(oldBackend) ){
				on.exit({
						if( verbose > 2 ){
							message("# Restoring previous foreach backend '", getDoBackendName(oldBackend) ,"' ... ", appendLF=FALSE)
						}
						setDoBackend(oldBackend, cleanup=TRUE)
						if( verbose > 2 ) message('OK')
					}, add=TRUE)
			}#
			
			# From this point, the backend is registered
			# => one knows if we'll run a sequential or parallel foreach loop
			.MODE_SEQ <- is.doSEQ()
			MODE_PAR <- .MODE_PAR <- !.MODE_SEQ
			
		}
		
		# check seed method: fixed values are not sensible -> warning
		.checkRandomness <- FALSE
		if( is.KINOMO(seed) && !is.empty.KINOMO(seed) ){
			.checkRandomness <- TRUE
		}
		# start_RNG_all
		# if the seed is numerical or a rstream object,	then use it to set the 
		# initial state of the random number generator:			
		# build a sequence of RNGstreams: if no suitable seed is provided
		# then the sequence use a random seed generated with a single draw 
		# of the current active RNG. If the seed is valid, then the 
		# 
		# setup the RNG sequence

		# override with standard RNG if .options$RNGstream=FALSE
		resetRNG <- NULL
		if( !.options$RNGstream && (!opt.parallel || .MODE_SEQ) ){
			
			.RNG.seed <- rep(list(NULL), nrun)
			if( isNumber(rng) ){
				resetRNG <- getRNG()
				if( verbose > 2 ) message("# Force using current RNG settings seeded with: ", rng)
				set.seed(rng)
			}else if( verbose > 2 ) 
				message("# Force using current RNG settings")
			
		}else{
			.RNG.seed <- setupRNG(rng, n = nrun, verbose=verbose)
			# restore the RNG state on exit as after RNGseq:
			# - if no seeding occured then the RNG has still been drawn once in RNGseq
			# which must be reflected so that different unseeded calls use different RNG states
			# - one needs to restore the RNG because it switched to L'Ecuyer-CMRG. 
			resetRNG <- getRNG()
		}
		stopifnot( length(.RNG.seed) == nrun )
		# update RNG settings on exit if necessary
		# and only if no error occured
		if( !is.null(resetRNG) ){
			on.exit({
				if( exitSuccess() ){
					if( verbose > 2 ) message("# Updating RNG settings ... ", appendLF=FALSE)							
					setRNG(resetRNG)
					if( verbose > 2 ) message("OK")
					if( verbose > 3 ) showRNG()
				}
			}, add=TRUE)
		}
		#end_RNG_all
		
		####FOREACH_KINOMO
		if( opt.parallel ){
			
			if( verbose ){
				if( verbose > 1 )
						message("# Using foreach backend: ", getDoParName()
								," [version ", getDoParVersion(),"]")
				# show number of processes
				if( getDoParWorkers() == 1 ) message("Mode: sequential [foreach:",getDoParName(),"]")
				else message("Mode: parallel ", str_c("(", getDoParWorkers(), '/', parallel::detectCores()," core(s))"))
			}
			
			# check shared memory capability
			.MODE_SHARED <- !keep.all && setupSharedMemory(verbose)
			
			# setup temporary directory when not keeping all fits
			if( !keep.all || verbose ){
				KINOMO_TMPDIR <- setupTempDirectory(verbose)
				# delete on exit
				if( .CLEANUP ){
					on.exit({
						if( verbose > 2 ) message("# Deleting temporary directory '", KINOMO_TMPDIR, "' ... ", appendLF=FALSE)
						unlink(KINOMO_TMPDIR, recursive=TRUE)
						if( verbose > 2 ) message('OK')
					}, add=TRUE)
				}
			}
			
			run.all <- function(x, rank, method, seed, model, .options, ...){
								
				## 1. SETUP
				# load some variables from parent environment to ensure they 
				# are exported in the foreach loop
				MODE_SEQ <- .MODE_SEQ
				MODE_SHARED <- .MODE_SHARED
				verbose <- verbose
				keep.all <- keep.all
				opt.gc <- .options$garbage.collect
				CALLBACK <- .callback
				.checkRandomness <- .checkRandomness
				
				# check if single or multiple host(s)
				hosts <- unique(getDoParHosts())
				if( verbose > 2 ) message("# Running on ", length(hosts), " host(s): ", str_out(hosts))
				SINGLE_HOST <- length(hosts) <= 1L
				MODE_SHARED <- MODE_SHARED && SINGLE_HOST
				if( verbose > 2 ) message("# Using shared memory ... ", MODE_SHARED)
				
				# setup mutex evaluation function
				mutex_eval <- if( MODE_SHARED ) ts_eval(verbose = verbose > 4) else force
				
				# Specific thing only if one wants only the best result
				if( !keep.all ){ 
					KINOMO_TMPDIR <- KINOMO_TMPDIR
					# - Define the shared memory objects
					vOBJECTIVE <- gVariable(as.numeric(NA), MODE_SHARED) 
					# the consensus matrix is computed only if not all the results are kept				
					vCONSENSUS <- gVariable(matrix(0, ncol(x), ncol(x)), MODE_SHARED)			
				}
				
				## 2. RUN
				# ensure that the package KINOMO is in each worker's search path
				.packages <- setupLibPaths('KINOMO', verbose>3)
                
                # export all packages that contribute to KINOMO registries, 
                # e.g., algorithms or seeding methods.
                # This is important so that these can be found in worker nodes
                # for non-fork clusters.
                if( !is.null(contribs <- registryContributors(package = 'KINOMO')) ){
                    .packages <- c(.packages, contribs)
                }
                
				# export dev environment if in dev mode 
#				.export <- if( isDevNamespace('KINOMO') && !is.doSEQ() ) ls(asNamespace('KINOMO'))
				
				# in parallel mode: verbose message from each run are only shown in debug mode
				.options$verbose <- FALSE 
				if( verbose ){
					if( debug || (.MODE_SEQ && verbose > 1) )
						.options$verbose <- verbose
					
					if( (!.MODE_SEQ && !debug) || (.MODE_SEQ && verbose == 1) ){
						if( verbose == 1 ){
							# create progress bar
							pbar <- txtProgressBar(0, nrun+1, width=50, style=3, title='Runs:'
													, shared=KINOMO_TMPDIR)
						}else{
							cat("Runs: ")
						}
					}
				}
				
				# get options from master process to pass to workers
				KINOMO.opts <- KINOMO.options()
				
				# load extra required packages for shared mode 
				if( MODE_SHARED ) 
					.packages <- c(.packages, 'bigmemory', 'synchronicity')
				
				res.runs <- foreach(n=1:nrun
								, RNGobj = .RNG.seed
								, .verbose = debug
								, .errorhandling = 'pass'
								, .packages = .packages
#								, .export = .export
#								, .options.RNG=.RNG.seed
								) %dopar% { #START_FOREACH_LOOP
				
					# Pass options from master process
					KINOMO.options(KINOMO.opts)
					
					# in mode sequential or debug: show details for each run
					if( MODE_SEQ && verbose > 1 )
						cat("\n## Run: ",n, "/", nrun, "\n", sep='')
					
					# set the RNG if necessary and restore after each run 
					if( MODE_SEQ && verbose > 2 )
							message("# Setting up loop RNG ... ", appendLF=FALSE)
					setRNG(RNGobj, verbose=verbose>3 && MODE_SEQ)
					if( MODE_SEQ && verbose > 2 )
							message("OK")
						
					# limited verbosity in simple mode
					if( verbose && !(MODE_SEQ && verbose > 1)){
						if( verbose >= 2 ) mutex_eval( cat('', n) )		
						else{
							# update progress bar (in mutex)
							mutex_eval(setTxtProgressBar(pbar, n))
							#
						}
					}
					
					# check RNG changes
					if( n == 1 && .checkRandomness ){
						.RNGinit <- getRNG()
					}
					
					# fit a single KINOMO model
					res <- KINOMO(x, rank, method, nrun=1, seed=seed, model=model, .options=.options, ...)
					
					if( n==1 && .checkRandomness && rng.equal(.RNGinit) ){
						warning("KINOMO::KINOMO - You are running multiple non-random KINOMO runs with a fixed seed")
					}
					
					# if only the best fit must be kept then update the shared objects
					if( !keep.all ){
						
						# initialise result list
						resList <- list(filename=NA, residuals=NA, .callback=NULL)
						
						##LOCK_MUTEX
						mutex_eval({
												
							# check if the run found a better fit
							.STATIC.err <- vOBJECTIVE()
							
							# retrieve approximation error
							err <- deviance(res)
							
							if( is.na(.STATIC.err) || err < .STATIC.err ){
								
								if( n>1 && verbose ){
									if( MODE_SEQ && verbose > 1 ) cat("## Better fit found [err=", err, "]\n")
									else if( verbose >= 2 ) cat('*')
								}
								
								# update residuals
								vOBJECTIVE(err)
								
								# update best fit on disk: use pid if not using shared memory
								resfile <- hostfile("fit", tmpdir=KINOMO_TMPDIR, fileext='.rds', pid=!MODE_SHARED)
								if( MODE_SEQ && verbose > 2 )
									message("# Serializing fit object in '", resfile, "' ... ", appendLF=FALSE)
								saveRDS(res, file=resfile, compress=FALSE)
								if( MODE_SEQ && verbose > 2 ){
									message(if( file.exists(resfile) ) 'OK' else 'ERROR')
								}
								# store the filename and achieved objective value in the result list
								resList$filename <- resfile
								resList$residuals <- err
																
							}
							
							## CONSENSUS
							# update the consensus matrix
							if( MODE_SHARED && SINGLE_HOST ){								
								# on single host: shared memory already contains consensus
								vCONSENSUS(vCONSENSUS() + connectivity(res, no.attrib=TRUE))
							}else{
								# on multiple hosts: must return connectivity and aggregate at the end
								resList$connectivity <- connectivity(res, no.attrib=TRUE)
							}
							
							## CALLBACK
							# call the callback function if necessary (return error as well)
							if( is.function(CALLBACK) ){
								resList$.callback <- tryCatch(CALLBACK(res, n), error=function(e) e)
							}
						
						})
						##END_LOCK_MUTEX
									
						# discard result object
						res <- NULL
						# return description list 
						res <- resList
					}
										
					# garbage collection if requested
					if( opt.gc && n %% opt.gc == 0 ){
						if( verbose > 2 ){
							if( MODE_SEQ )
								message("# Call garbage collector")
							else{
								mutex_eval( cat('%') )
							}
						}
						
						gc(verbose= MODE_SEQ && verbose > 3)
					}
					
					# return the result
					res
				}				
				## END_FOREACH_LOOP
				
				if( verbose && !debug ){
					if( verbose >= 2 ) cat(" ... DONE\n")
					else{
						setTxtProgressBar(pbar, nrun+1)
						pbar$kill(.CLEANUP)
					}
				}
				
				## 3. CHECK FIT ERRORS
				errors <- checkErrors(res.runs)
				if( errors$n > 0 ){
					fstop(errors$n,"/", nrun, " fit(s) threw an error.\n"
							,"# Error(s) thrown:\n", errors$msg)
				}
				
				## 4. WRAP UP
				if( keep.all ){ # result is a list of fits
					# directly return the list of fits
					res <- res.runs
					
				}else{ # result is a list of lists: filename, .callback 
					# loop over the result files to find the best fit
					if( verbose > 2 ) message("# Processing partial results ... ", appendLF=FALSE)
					ffstop <- function(...){ message('ERROR'); fstop(...) }
					# get best fit index
					idx <- which.min(sapply(res.runs, '[[', 'residuals'))
					if( length(idx) == 0L )
						ffstop("Unexpected error: no partial result seem to have been saved.")
					resfile <- res.runs[[idx]]$filename
					# check existence of the result file
					if( !file_test('-f', resfile) )
						ffstop("could not find temporary result file '", resfile, "'")
						
					# update res with a better fit
					res <- readRDS(resfile)
					if( !isKINOMOfit(res) ) 
						ffstop("invalid object found in result file '", resfile, "'")
					if( verbose > 2 ) message('OK')
					# wrap the result in a list: fit + consensus
					res <- list(fit=res, consensus=NA)
					
					# CONSENSUS MATRIX
					if( !is.null(res.runs[[1]]$connectivity) ){ # not MODE_SHARED
						# aggregate connectivity matrices
						con <- matrix(0, ncol(x), ncol(x))
						sapply(res.runs, function(x){
							con <<- con + x$connectivity 
						})
						res$consensus <- con
						
					}else{ # in MODE_SHARED: get consensus from global shared variable
						res$consensus <- vCONSENSUS()
						cn <- colnames(x) 
						if( is.null(cn) ) dimnames(res$consensus) <- NULL
						else dimnames(res$consensus) <- list(cn, cn)
					}
					
					# CALLBACKS
					if( !is.null(.callback) ){
						res$.callback <- processCallback(res.runs)
					}
				}
				##
				
				if( MODE_SEQ && verbose>1 ) cat("## DONE\n")
				
				# return result
				res
			}			
		}####END_FOREACH_KINOMO
		else{####SAPPLY_KINOMO
			
			run.all <- function(x, rank, method, seed, model, .options, ...){
				
				# by default force no verbosity from the runs
				.options$verbose <- FALSE
				if( verbose ){
					message("Mode: sequential [sapply]")
					if( verbose > 1 ){					
						# pass verbosity options in this case
						.options$verbose <- verbose
					}
				}
				
				## 1. SETUP				
				# define static variables for the case one only wants the best result
				if( !keep.all ){
					# statis list with best result: fit, residual, consensus
					best.static <- list(fit=NULL, residuals=NA, consensus=matrix(0, ncol(x), ncol(x)))					
				}
											
				## 2. RUN:
				# perform a single run `nrun` times
				if( verbose == 2 ){
					showRNG()
				}
				if( verbose && !debug ) cat('Runs:')
				res.runs <- mapply(1:nrun, .RNG.seed, FUN=function(n, RNGobj){
					
					#start_verbose
					if( verbose ){
						# in mode verbose > 1: show details for each run
						if( verbose > 1 ){
							cat("\n## Run: ",n, "/", nrun, "\n", sep='')							
						}else{
						# otherwise only some details for the first run
							cat('', n)
						}
					}#end_verbose
					
					# set the RNG for each run
					if( verbose > 2 ) message("# Setting up loop RNG ... ", appendLF=FALSE)
					setRNG(RNGobj, verbose=verbose>3)
					if( verbose > 2 ) message("OK")
					
					# check RNG changes
					if( n == 1 && .checkRandomness ){
						.RNGinit <- getRNG()
					}
					
					# fit a single KINOMO model
					res <- KINOMO(x, rank, method, nrun=1, seed=seed, model=model, .options=.options, ...)
					
					if( n==1 && .checkRandomness && rng.equal(.RNGinit) ){
						warning("KINOMO::KINOMO - You are running multiple non-random KINOMO runs with a fixed seed"
								, immediate.=TRUE)
					}
					
					if( !keep.all ){
						
						# initialise result list
						resList <- list(residuals=NA, .callback=NULL)
						
						# check if the run found a better fit
						err <- residuals(res)
						best <- best.static$residuals
						if( is.na(best) || err < best ){
							if( verbose ){
								if( verbose > 1L ) cat("## Updating best fit [deviance =", err, "]\n", sep='')
								else cat('*')
							}
							
							# update best fit (only if necessary)
							best.static$fit <<- res
							best.static$residuals <<- err
							
							resList$residuals <- err
						}
							
						# update the static consensus matrix (only if necessary)
						best.static$consensus <<- best.static$consensus + connectivity(res, no.attrib=TRUE)
						
						# call the callback function if necessary
						if( !is.null(.callback) ){
							resList$.callback <- tryCatch(.callback(res, n), error=function(e) e)							
						}
						
						# reset the result to NULL
						res <- resList
						
					}
					
					# garbage collection if requested
					if( opt.gc && n %% opt.gc == 0 ){
						if( verbose > 1 )
							message("# Call garbage collection NOW")
						else if( verbose )
							cat('%')
						
						gc(verbose = verbose > 3)
					}
					
					if( verbose > 1 ) cat("## DONE\n")
					
					# return the result
					res
				}, SIMPLIFY=FALSE)
				##
				
				if( verbose && !debug ) cat(" ... DONE\n")
				
				## 3. ERROR CHECK / WRAP UP
				
				if( keep.all ){
					res <- res.runs
				}else{
					res <- list(fit=best.static$fit, consensus=best.static$consensus)
					
					# CALLBACKS
					if( !is.null(.callback) ){
						res$.callback <- processCallback(res.runs)
					}
				}
				
				res
			}
			
		}####END_SAPPLY_KINOMO
			
		####END_DEFINE_RUN		
		
		# perform all the KINOMO runs
		t <- system.time({res <- run.all(x=x, rank=rank, method=method, seed=seed, model=model, .options, ...)})
		if( verbose && !debug ){
			cat("System time:\n")
			print(t)
		}
		
		if( keep.all ){
			
			# when keeping all the fits: join the results into an KINOMOfitXn object
			# TODO: improve memory management here
			res <- KINOMOfitX(res, runtime.all=t)
			
			return( exitSuccess(res) )
			
		}else{# if one just want the best result only return the best
			# ASSERT the presence of the result
			stopifnot( !is.null(res$fit) )
			# ASSERT the presence of the consensus matrix
			stopifnot( !is.null(res$consensus) )
			
			res.final <- KINOMOfitX(res$fit, consensus=res$consensus/nrun
								, runtime.all=t, nrun=as.integer(nrun)
								, rng1=.RNG.seed[[1]])
			
			# ASSERT and add callback if necessary 
			if( !is.null(.callback) ){
				stopifnot( !is.null(res$.callback) )
				res.final$.callback <- res$.callback
			}
			
			return( exitSuccess(res.final) )
		}
		
	}##END_MULTI_RUN
	
	# start_RNG
	# show original RNG settings in verbose > 2
	if( verbose > 3 ){
		message("# ** Current RNG settings:")
		showRNG()
	}
	
	# do something if the RNG was actually changed
	newRNG <- getRNG()
	.RNG.seed <- setupRNG(rng, 1, verbose=verbose-1)
	# setup restoration
	if( isRNGseed(rng) ){
		if( verbose > 3 ) showRNG()
		
		# restore RNG settings
		on.exit({
			if( verbose > 2 ) message("# Restoring RNG settings ... ", appendLF=FALSE)							
			setRNG(newRNG)					
			if( verbose > 2 ) message("OK")
			if( verbose > 3 ) showRNG()
		}, add=TRUE)
	}
	#end_RNG
	
	# CHECK PARAMETERS:	
	# test for negative values in x only if the method is not mixed
	if( !is.mixed(method) && min(x, na.rm = TRUE) < 0 )
        fstop('Input matrix ', substitute(x),' contains some negative entries.');
	# test if one row contains only zero entries
    if( min(rowSums(x, na.rm = TRUE), na.rm = TRUE) == 0 )
        fstop('Input matrix ', substitute(x),' contains at least one null or NA-filled row.');	

	# a priori the parameters for the run are all the one in '...'
	# => expand with the strategy's defaults (e.g., maxIter)
	parameters.method <- expand_list(list(...), .method_defaults)
	#
	
	if( is.KINOMO(seed) ){
		
		if( !is.null(model) )
			fwarning("Discarding argument `model`: directly using KINOMO model supplied in argument `seed`")
		
		# if the seed is a KINOMOfit object then only use the fit (i.e. the KINOMO model)
		# => we want a fresh and clean KINOMOfit object
		if( isKINOMOfit(seed) )
			seed <- fit(seed)
		
		# Wrap up the seed into a KINOMOfit object
		seed <- KINOMOfit(fit=seed, seed='KINOMO')
	}
	else if( !inherits(seed, 'KINOMOfit') ){
		
		## MODEL INSTANTIATION :
	
		# default KINOMO model is retrieved from the KINOMO strategy
		.modelClass <- modelname(method)
		# if a character string then use this type of KINOMO model, but still look 
		# for slots in `...`
		if( is.character(model) ){
			.modelClass <- model
			model <- NULL
		}
		
		# some of the instantiation parameters are set internally
		# TODO: change target into x (=> impact on KINOMOModel ?
		parameters.model.internal <- list(rank=rank, target=0)
		parameters.model <- list()
		
		init <- 
		if( is.KINOMO(model) ){
			model
		}else{
			# if 'model' is NULL: initialization parameters are searched in '...' 
			if( is.null(model) ){
				
				# extract the parameters from '...' that correspond to slots in the given class
				stopifnot( isKINOMOclass(.modelClass) )
				parameters <- .extract.slots.parameters(.modelClass, parameters.method)	
				
				# restrict parameters.method to the ones that won't be used to instantiate the model
				overriden <- is.element(names(parameters$slots), names(parameters.model.internal))
				parameters.method <- c(parameters$extra, parameters$slots[overriden])
							
				#- the model parameters come from the remaining elements
				parameters.model <- c(model=.modelClass, parameters$slots)
				
			} else if( is.list(model) ){  # otherwise argument 'model' must be a list
				
				# if the list is not empty then check all elements are named and 
				# not conflicting with the internally set values			
				if( length(model) > 0 ){
					# all the elements must be named
					if( !hasNames(model, all=TRUE) )  
						fstop("Invalid argument `model` [elements must all be named]. See ?KINOMO.")
					
					# warn the user if some elements are conflicting and won't be used
					overriden <- is.element(names(model), names(parameters.model.internal))
					if( any(overriden) )
						warning("KINOMO::KINOMO - Model parameter(s) [" 
								, str_out(model[overriden], use.names=TRUE, max=Inf)
								, "] discarded. Used internally set value(s) ["
								, str_out(parameters.model.internal[names(model[overriden])], use.names=TRUE, max=Inf)
								, "]"
								, call.=FALSE)
				}
				
				# add default model class if necessary
				if( is.null(model$model) )
					model$model <- .modelClass
				# all the instantiation parameters come from argument 'model'
				parameters.model <- model
				
			}else{ 			
				fstop("Invalid argument 'model' [expected NULL, a character string, or a list to set slots in the KINOMO model class '",.modelClass,"']. See ?KINOMO.")
			}	
				
			
			#- force the value of the internally set arguments for the instantiation of the model
			parameters.model <- .merge.override(parameters.model, parameters.model.internal)		
			
			# at this point 'init' should be the list of the initialization parameters
			if( !is.list(parameters.model) ){
				fstop("Unexpected error: object 'parameters.model' must be a list")
			}
			if( !is.element('model', names(parameters.model)) ){
				fstop("Unexpected error: object 'parameters.model' must contain an element named 'model'")
			}
			
			parameters.model
		}
	## SEEDING:
	# the seed must either be an instance of class 'KINOMO', the name of a seeding method as a character string
	# or a list of parameters to pass to the 'seed' function.
			parameters.seed <- list()
			seed.method <- NULL
			if( (is.character(seed) && length(seed) == 1) 
				|| is.numeric(seed) 
				|| is.null(seed) 
#				|| is(seed, 'rstream') 
				) seed.method <- seed
			else if( is.function(seed) ) seed.method <- seed
			else if( is.list(seed) ){ # seed is a list...
				
				if( !is.null(seed$method) ){ # 'seed' must contain an element giving the method...
					seed.method <- seed$method
					parameters.seed <- seed[-which(names(seed)=='method')]
				}
				else if ( is.null(names(seed)) || names(seed)[1] == '' ){ # ... or the first element must be a method
					seed.method <- seed[[1]]
					if( length(seed) > 1 ) parameters.seed <- seed[2:length(seed)]
				}
				else fstop("Invalid parameter: list 'seed' must contain the seeding method through its first element or through an element named 'method' [", str_desc(seed, 2L), "]")
				
				# check validity of the method provided via the list
				if( !is.function(seed.method) && !(is.character(seed.method) && length(seed.method)==1) )
					fstop("The seeding method provided by parameter 'seed' [", str_desc(seed.method), "] is invalid: a valid function or a character string is expected")
			}
			else fstop("Invalid parameter 'seed'. Acceptable values are:\n\t- ",
						paste("an object that inherits from class 'KINOMO'"
							, "the name of a seeding method (see ?KINOMOSeed)"
							, "a valid seed method definition"
							, "a list containing the seeding method (i.e. a function or a character string) as its first element\n\tor as an element named 'method' [and optionnally extra arguments it will be called with]"
							, "a numerical value used to set the seed of the random generator"
							, "NULL to directly pass the model instanciated from arguments 'model' or '...'."
							, sep="\n\t- "))
						 			
			# call the 'seed' function passing the necessary parameters
			if( verbose )
				message("KINOMO seeding method: ", 
						if( is.character(seed.method) || is.numeric(seed.method) ) seed.method
						else if( is.null(seed.method) ) 'NULL'
						else if( !is.null(attr(seed.method, 'name')) ) attr(seed.method, 'name') 
						else if( is.function(seed.method) ) '<function>'
						else NA)
			
			#seed <- do.call(getGeneric('seed', package='KINOMO')
			seed <- do.call(getGeneric('seed')
					, c(list(x=x, model=init, method=seed.method), parameters.seed))
			
			# check the validity of the seed
			if( !inherits(seed, 'KINOMOfit') ) 
				fstop("The seeding method function should return class 'KINOMO' ["
					, if( is.character(seed.method) ) paste('method "', seed.method, "' ", sep='') else NULL 
					, "returned class: '", class(seed), "']")
	}
	# -> at this point the 'seed' object is an instance of class 'KINOMOfit'
	KINOMO.debug('KINOMO', "Seed is of class: '", class(seed), "'")
	# ASSERT just to be sure
	if( !inherits(seed, 'KINOMOfit') )
		fstop("Invalid class '", class(seed), "' for the computed seed: object that inherits from class 'KINOMOfit' expected.")
	
	# check the consistency of the KINOMO model expected by the algorithm and 
	# the one defined by the seed
	#if( none( sapply(model(method), function(c) extends(model(seed), c)) ) )
	if( all( !inherits(fit(seed), modelname(method)) ) )
		fstop("Invalid KINOMO model '", modelname(seed),"': algorithm '", name(method), "' expects model(s) "
			, paste(paste("'", modelname(method),"'", sep=''), collapse=', ')
			, " or extension.")
	
	# get the complete seeding method's name 
	seed.method <- seeding(seed)
	
	## FINISH SETUP OF THE SEED OBJECT: store some data within the seed so
	# that strategy methods can access them directly
	algorithm(seed) <- name(method) # algorithm name
	seed@distance <- objective(method) # distance name
	seed@parameters <- parameters.method # extra parameters
	run.options(seed) <- KINOMO.options() # set default run options
	run.options(seed, 'error.track') <- .OPTIONS$track
	if( is.numeric(.OPTIONS$track) )
		run.options(seed, 'track.interval') <- .OPTIONS$track
	run.options(seed, 'verbose') <- verbose
	# store ultimate KINOMO() call
	seed@call <- match.call()
	##

	## print options if in verbose > 3
	if( verbose > 3 ){
		cat("## OPTIONS:\n")		
		sapply(seq_along(.options)
				, function(i){
					r <- i %% 4
					cat(if(r!=1) '\t| ' else "# ", names(.options)[i],': ', .options[[i]], sep='')
					if(r==0) cat("\n")
				})
		if( length(.options) %% 4 != 0 )cat("\n")
	}
	
	
	## run parameters: 
	parameters.run <- c(list(object=method, y=x, x=seed), parameters.method)
	## Compute the initial residuals if tracking is enabled
	init.resid <- if( .OPTIONS$track && !is.partial.KINOMO(seed) ){
		do.call('deviance', parameters.run)
	}
	
	## RUN KINOMO METHOD:
	# call the strategy's run method [and time it]
	t <- system.time({				
		res <- if( !dry.run ){
				do.call('run', parameters.run)
				
			}else{ 
				seed
			}
	})

	## WRAP/CHECK RESULT
	res <- .wrapResult(x, res, seed, method=method, seed.method=seed.method, t)
	if( !isKINOMOfit(res) ){ # stop if error
		fstop(res)
	}
	##
	
	## CLEAN-UP + EXTRAS:
	# add extra information to the object
	# slot 'parameters'
	if( length(res@parameters) == 0L && length(parameters.method)>0L )
		res@parameters <- parameters.method
	# last residuals
	if( length(residuals(res)) == 0 && !is.partial.KINOMO(seed) ){
		parameters.run$x <- res
		residuals(res, niter=niter(res)) <- do.call('deviance', parameters.run)
	}
	# first residual if tracking is enabled
	if( .OPTIONS$track && !is.null(init.resid) ){
		if( !hasTrack(res, niter=0) )
			residuals(res, track=TRUE) <- c('0'=init.resid, residuals(res, track=TRUE))
	}
	
	if( length(residuals(res)) && is.na(residuals(res)) ) warning("KINOMO residuals: final objective value is NA")
	res@runtime <- t
	
	# return the result
	exitSuccess(res)
})

# wrap result
.wrapResult <- function(x, res, seed, method, seed.method, t){
	
	## wrap into an KINOMOfit object (update seed)
	if( !isKINOMOfit(res) ){
		# extract expression data if necessary
		if( is(res, 'ExpressionSet') ) res <- exprs(res)
		if( is(x, 'ExpressionSet') ) x <- exprs(x)
		# wrap
		if( is.matrix(res) ){
			if( ncol(res) == ncol(x) ){# partial fit: coef
				# force dimnames
				colnames(res) <- colnames(x)
				res <- KINOMOModel(H=res)
			}else if( nrow(res) == nrow(x) ){# partial fit: basis
				# force dimnames
				rownames(res) <- rownames(x)
				res <- KINOMOModel(W=res)
			}
		}else if( is.list(res) ){ # build KINOMO model from result list
			res <- do.call('KINOMOModel', res)
		}
		
		# substitute model in fit object
		if( is.KINOMO(res) ){
			tmp <- seed
			fit(tmp) <- res
			tmp@runtime <- t
			res <- tmp
		}
	}
	
	## check result
	if( !isTRUE(err <- .checkResult(res, seed)) ) return(err) 
	
	## Enforce some slot values
	# slot 'method'
	algorithm(res) <- name(method)	
	# slot 'distance'
	res@distance <- objective(method)	
	# slot 'seed'
	if( seed.method != '' ) seeding(res) <- seed.method
	# set dimnames of the result only if necessary
	if( is.null(dimnames(res)) )
		dimnames(res) <- dimnames(seed)
	
	res
}

# check result
.checkResult <- function(fit, seed){
	# check the result is of the right type
	if( !inherits(fit, 'KINOMOfit') ){ 
		return(str_c("KINOMO algorithms should return an instance of class 'KINOMOfit' [returned class:", class(fit), "]"))
	}
	
	# check that the model has been fully estimated
	if( is.partial.KINOMO(fit) ){
		warning("KINOMO - The KINOMO model was only partially estimated [dim = (", str_out(dim(fit), Inf),")].")
	}
	# check that the fit conserved all fixed terms (only warning)
	if( nterms(seed) ){
		if( length(i <- icterms(seed)) && !identical(coef(fit)[i,], coef(seed)[i,]) ){
			warning("KINOMO - Fixed coefficient terms were not all conserved in the fit: the method might not support them.")
		}
		if( length(i <- ibterms(seed)) && !identical(basis(fit)[,i], basis(seed)[,i]) ){
			warning("KINOMO - Fixed basis terms were not all conserved in the fit: the method might not support them.")
		}
	}
	TRUE
}


setGeneric('seed', function(x, model, method, ...) standardGeneric('seed') )


setMethod('seed', signature(x='matrix', model='KINOMO', method='KINOMOSeed'), 
	function(x, model, method, rng, ...){	
		
		# debug message
		KINOMO.debug('seed', "use seeding method: '", name(method), "'")
		
		# temporarly set the RNG if provided 
		if( !missing(rng) ){
			orng <- setRNG(rng)
			on.exit(setRNG(orng))
		}
		
		# save the current RNG numerical seed
		rng.s <- getRNG()
		# create the result KINOMOfit object, storing the RNG numerical seed
		res <- KINOMOfit()
		# ASSERT: check that the RNG seed is correctly set
		stopifnot( rng.equal(res,rng.s) )
		# call the seeding function passing the extra parameters
		f <- do.call(algorithm(method), c(list(model, x), ...))
		# set the dimnames from the target matrix
		dimnames(f) <- dimnames(x)
		# set the basis names from the model if any
		if( !is.null(basisnames(model)) )
			basisnames(f) <- basisnames(model)
		# store the result into the KINOMOfit object
		fit(res) <- f
		
		# if not already set: store the seeding method's name in the resulting object
		if( seeding(res) == '' ) seeding(res) <- name(method)
		
		# return the seeded object
		res
	}
)

setMethod('seed', signature(x='ANY', model='ANY', method='function'),
	function(x, model, method, name, ...){
		
		# generate runtime name if necessary
		if( missing(name) ) name <- basename(tempfile("KINOMO.seed."))
		# check that the name is not a registered name		
		if( existsKINOMOSeed(name) )
			stop("Invalid name for custom seeding method: '",name,"' is already a registered seeding method")
		
		# wrap function method into a new KINOMOSeed object		 						
		seedObj <- new('KINOMOSeed', name=name, method=method)
		# call version with KINOMOSeed 
		seed(x, model, seedObj, ...)
	}
)


setMethod('seed', signature(x='ANY', model='ANY', method='missing'),
	function(x, model, method, ...){
		seed(x, model, KINOMO.getOption('default.seed'), ...)
	}
)

setMethod('seed', signature(x='ANY', model='ANY', method='NULL'),
	function(x, model, method, ...){
		seed(x, model, 'none', ...)
	}
)

setMethod('seed', signature(x='ANY', model='ANY', method='numeric'),
	function(x, model, method, ...){
		
		# set the seed using the numerical value by argument 'method'
		orng <- setRNG(method)
		#TODO: restore the RNG state?
		
		# call seeding method 'random'
		res <- seed(x, model, 'random', ...)
		
		# return result
		return(res)
	}
)

setMethod('seed', signature(x='ANY', model='ANY', method='character'),
		function(x, model, method, ...){
			
			# get the seeding method from the registry
			seeding.fun <- KINOMOSeed(method)
				
			#Vc#Use seeding method: '${method}'
			# call 'seed' with the seeding.function			
			seed(x, model, method=seeding.fun, ...)
			
		}
)

setMethod('seed', signature(x='ANY', model='list', method='KINOMOSeed'), 
	function(x, model, method, ...){	
		
		## check validity of the list: there should be at least the KINOMO (sub)class name and the rank
		if( length(model) < 2 )
			stop("Invalid parameter: list 'model' must contain at least two elements giving the model's class name and the factorization rank")
					
		# 'model' must contain an element giving the class to instanciate
		if( is.null(model$model) ){
			
			err.msg <- "Invalid parameter: list 'model' must contain a valid KINOMO model classname in an element named 'model' or in its first un-named element"			
			unamed <- if( !is.null(names(model)) ) which(names(model) %in% c('', NA)) else 1			
			if ( length(unamed) > 0 ){ # if not the first unamed element is taken as the class name
				idx <- unamed[1]
				val <- unlist(model[idx], recursive=FALSE)				
				if( is.character(val) && length(val)==1 && extends(val, 'KINOMO') )
					names(model)[idx] <- 'model'
				else stop(err.msg)
			}else stop(err.msg)
		}
		
		# 'model' must contain an element giving the factorization rank
		if( is.null(model$rank) ){
			err.msg <- "Invalid parameter: list 'model' must contain the factorization rank in an element named 'rank' or in its second un-named element"
			unamed <- if( !is.null(names(model)) ) which(names(model) %in% c('', NA)) else 1
			if ( length(unamed) > 0 ){ # if not the second element is taken as the factorization rank
				idx <- unamed[1]
				val <- unlist(model[idx], recursive=FALSE)
				if( is.numeric(val) && length(val)==1 )
					names(model)[idx] <- 'rank'
				else stop(err.msg)
			}
			else stop(err.msg)
		}
					
		KINOMO.debug('seed', "using model parameters:\n", capture.output(print(model)) )
		# instantiate the object using the factory method		
		model <- do.call('KINOMOModel', model)
		KINOMO.debug('seed', "using KINOMO model '", class(model), "'")
		
		# check that model is from the right type, i.e. inherits from class KINOMO
		if( !inherits(model, 'KINOMO') ) stop("Invalid object returned by model: object must inherit from class 'KINOMO'")
		
		seed(x, model, method, ...)
	}
)

setMethod('seed', signature(x='ANY', model='numeric', method='KINOMOSeed'), 
	function(x, model, method, ...){	

		seed(x, KINOMOModel(model), method, ...)
	}
)



.extract.slots.parameters <- function(class.name, ...){
		
	# check validity of class.name
	if( !isClass(class.name) ) stop("Invalid class name: class '", class.name, "' dose not exist")
	
	# transform '...' into a list
	parameters <- list(...)	
	
	if( length(parameters) == 1L && is.null(names(parameters)) ){
		parameters <- parameters[[1L]]
	}

	# get the slots from the class name
	slots <- slotNames(class.name)
	# get the named parameters that correspond to a slot
	in.slots <- is.element(names(parameters), slots)
	# return the two lists	
	list( slots=parameters[in.slots], extra=parameters[!in.slots])
}


.merge.override <- function(l1, l2, warning=FALSE){
	sapply(names(l2), function(name){
				if( warning && !is.null(l1[[name]]) )
					warning("overriding element '", name, "'")
				l1[[name]] <<- l2[[name]]
			})
	
	# return updated list
	return(l1)
}


KINOMOEstimateRank <- function(x, range, method=KINOMO.getOption('default.algorithm')
					, nrun=30, model=NULL, ..., verbose=FALSE, stop=FALSE){
	
	# fix method if passed NULL (e.g., from KINOMO('formula', 'numeric'))
	if( is.null(method) )
		method <- KINOMO.getOption('default.algorithm')
	
	# special handling of formula: get target data from the formula 
	if( is(x, 'formula') ){
		# dummy model to resolve formula
		dummy <- KINOMOModel(x, 0L, data=model)
		# retrieve target data
		V <- attr(dummy, 'target')
	}else{
		V <- x
	}
	
	# remove duplicates and sort
	range <- sort(unique(range))
	
	# initiate the list of consensus matrices: start with single NA values
	c.matrices <- setNames(lapply(range, function(x) NA), as.character(range))
	fit <- setNames(lapply(range, function(x) NA), as.character(range))
	bootstrap.measures <- list()

	# combine function: take all the results at once and merge them into a big matrix
	comb <- function(...){
		measures <- list(...)
		
		err <- which( sapply(measures, is.character) )		
		if( length(err) == length(measures) ){ # all runs produced an error
		
			# build an warning using the error messages
			msg <- paste(paste("#", seq_along(range),' ', measures, sep=''), collapse="\n\t-")
			stop("All the runs produced an error:\n\t-", msg)
		
		}else if( length(err) > 0 ){ # some of the runs returned an error
			
			# simplify the results with no errors into a matrix
			measures.ok <- sapply(measures[-err], function(x) x)
			
			# build a NA matrix for all the results
			n <- nrow(measures.ok)
			tmp.res <- matrix(as.numeric(NA), n, length(range))
			rownames(tmp.res) <- rownames(measures.ok)
			
			# set the results that are ok
			tmp.res[,-err] <- measures.ok
			# set only the rank for the error results 
			tmp.res['rank', err] <- range[err]
			# build an warning using the error messages
			msg <- paste(paste("#", err, measures[err], ' ', sep=''), collapse="\n\t-")
			warning("NAs were produced due to errors in some of the runs:\n\t-", msg)
			
			# return full matrix
			tmp.res
		}
		else # all the runs are ok 
			sapply(measures, function(x) x)
	}
	
#	measures <- foreach(r = range, .combine=comb, .multicombine=TRUE, .errorhandling='stop') %do% {
	k.rank <- 0 
	measures <- sapply(range, function(r, ...){
			k.rank <<- k.rank + 1L
			if( verbose ) cat("Compute KINOMO rank=", r, " ... ")
			
			# restore RNG on exit (except after last rank)
			# => this ensures the methods use the same stochastic environment
			orng <- RNGseed()
			if( k.rank < length(range) ) on.exit( RNGseed(orng), add = TRUE)
			
			res <- tryCatch({ #START_TRY
				
				res <- KINOMO(x, r, method, nrun=nrun, model=model, ...)					
				# directly return the result if a valid KINOMO result
				if( !isKINOMOfit(res, recursive = FALSE) )
					return(res)
				
				# store the consensus matrix
				c.matrices[[as.character(r)]] <<- consensus(res)
				# store the fit
				fit[[as.character(r)]] <<- res				
				
				# if confidence intervals must be computed then do it
	#			if( conf.interval ){
	#				# resample the tries
	#				samp <- sapply(seq(5*nrun), function(i){ sample(nrun, nrun, replace=TRUE) })
	#				
	#				bootstrap.measures[[as.character(r)]] <<- apply(samp, 2, function(s){
	#					res.sample <- join(res[s])
	#					summary(res.sample, target=x)
	#				})
	#			}
				
				# compute quality measures
				if( verbose ) cat('+ measures ... ')
				measures <- summary(res, target=V)
				
				if( verbose ) cat("OK\n")
				
				# return the measures
				measures
			} #END_TRY
	
			, error = function(e) {
					mess <- if( is.null(e$call) ) e$message else paste(e$message, " [in call to '", e$call[1],"']", sep='')
					mess <- paste('[r=', r, '] -> ', mess, sep='')
					if( stop ){ # throw the error
						if( verbose ) cat("\n")
						stop(mess, call.=FALSE)
					} # pass the error message
					if( verbose ) message("ERROR")					
					return(mess)
				}
			)
			
			# return the result
			res
		}
	, ..., simplify=FALSE)
	
	measures <- do.call(comb, measures)
	
	# reformat the result into a data.frame
	measures <- as.data.frame(t(measures))	
	
	# wrap-up result into a 'KINOMO.rank' S3 object
	res <- list(measures=measures, consensus=c.matrices, fit=fit)
	#if( conf.interval ) res$bootstrap.measure <- bootstrap.measures
	class(res) <- 'KINOMO.rank'
	return(res)
	
}


summary.KINOMO.rank <- function(object, ...){
	s <- summary(new('KINOMOList', object$fit), ...)
	# NB: sort measures in the same order as required in ...
	i <- which(!names(s) %in% names(object$measures))
	cbind(s[, i], object$measures[match(object$measures$rank, s$rank), ])
}



plot.KINOMO.rank <- function(x, y=NULL, what=c('all', 'cophenetic', 'rss', 'residuals'
									, 'dispersion', 'evar', 'sparseness'
									, 'sparseness.basis', 'sparseness.coef'
                                    , 'silhouette'
                                    , 'silhouette.coef', 'silhouette.basis'
                                    , 'silhouette.consensus')
						, na.rm=FALSE
                        , xname = 'x'
                        , yname = 'y'
                        , xlab = 'Factorization rank'
                        , ylab = ''
                        , main = 'KINOMO rank survey'
                        , ... ){

	
    # trick for convenience 
	if( is.character(y) && missing(what) ){
		what <- y
		y <- NULL
	}
	
	what <- match.arg(what, several.ok=TRUE)
    if( 'all' %in% what ){
        what <- c('cophenetic', 'rss', 'residuals', 'dispersion', 'evar', 'sparseness', 'silhouette')
    }
    
    .getvals <- function(x, xname){
    	measures <- x$measures
    	iwhat <- unlist(lapply(paste('^',what,sep=''), grep, colnames(measures)))
    	
    	# remove NA values if required
    	if( na.rm )
    		measures <- measures[ apply(measures, 1, function(row) !any(is.na(row[iwhat]))), ]
    	
    	vals <- measures[,iwhat, drop=FALSE]
    	x <- as.numeric(measures$rank)
    	xlim <- range(x)
        
        # define measure type
        measure.type <- setNames(rep('Best fit', ncol(measures)), colnames(measures))
        cons.measures <- c('silhouette.consensus', 'cophenetic', 'cpu.all')
        measure.type[match(cons.measures, names(measure.type))] <- 'Consensus'
        measure.type[grep("\\.coef$", names(measure.type))] <- 'Coefficients'
        measure.type[grep("\\.basis$", names(measure.type))] <- 'Basis'
        measure.type <- factor(measure.type)
        
        pdata <- melt(cbind(rank = x, vals), id.vars = 'rank')
        # set measure type
        pdata$Type <- measure.type[as.character(pdata$variable)]
        # define measure groups
        pdata$Measure <- gsub("^([^.]+).*", "\\1", pdata$variable)
        pdata$Data <- xname
        pdata
    }
    
    pdata <- .getvals(x, xname)
    
    # add reference data
    if( is(y, 'KINOMO.rank') ){
        pdata.y <- .getvals(y, yname)
        pdata <- rbind(pdata, pdata.y)
    }
    
    p <- ggplot(pdata, aes_string(x = 'rank', y = 'value')) +
            geom_line( aes_string(linetype = 'Data', colour = 'Type') ) +
            geom_point(size = 2, aes_string(shape = 'Data', colour = 'Type') ) +
            theme_bw() +
            scale_x_continuous(xlab, breaks = unique(pdata$rank)) +
            scale_y_continuous(ylab) +
            ggtitle(main)
    # remove legend if not necessary
    if( !is(y, 'KINOMO.rank') ){
        p <- p + scale_shape(guide = 'none') + scale_linetype(guide = 'none')
    }
    
    # use fix set of colors
    myColors <- brewer.pal(5,"Set1")
    names(myColors) <- levels(pdata$Type)
    p <- p + scale_colour_manual(name = "Measure type", values = myColors)
    
    # add facet
    p <- p + facet_wrap( ~ Measure, scales = 'free')
    
    # return plot
    p
}

