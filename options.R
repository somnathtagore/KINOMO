

#' @name options-KINOMO
NULL
.OPTIONS <- setupPackageOptions(
	# default algorithm
	default.algorithm='brunet'
	# default seeding method
	, default.seed='random'
	# track error during KINOMO updates
	, error.track = option_symlink('track') # for backward compatibility
	, track=FALSE
	# define the tracking interval
	, track.interval=30
	# define garbage collection interval
	, gc=50
	# define default parallel backend 
	, parallel.backend= option_symlink('pbackend') # for backward compatibility
	, pbackend= if( parallel::detectCores() > 1 ) 'par' else 'seq'
	# toogle verbosity
	, verbose=FALSE
	# toogle debug mode
	, debug=FALSE
, RESET=TRUE)


#' 
KINOMO.options <- .OPTIONS$options


KINOMO.getOption <- .OPTIONS$getOption


KINOMO.resetOptions <- .OPTIONS$resetOptions


KINOMO.printOptions <- .OPTIONS$printOptions



# debugging utility
KINOMO.debug <- function(fun, ...){
	if( KINOMO.getOption('debug') ){
		call.stack <- sys.calls()
		n <- length(call.stack)
		if( is.null(fun) ) fun <- as.character(call.stack[[n-1]]) 
		message('DEBUG::', fun, ' -> ', ...)
	}
	return(invisible())
}
