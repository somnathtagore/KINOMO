


#' @include algorithmic.R
#' @include KINOMOSet-class.R
NULL


setClass('Strategy'
	, contains = 'VIRTUAL'
	, representation = representation(
		name = 'character' # the strategy name
		, package = 'character' # the package that defines the strategy
		, defaults = 'list'
	)
	, prototype = prototype(
		package = character()
		, name = character()
	)
	, validity=function(object){
		
		# slot 'name' must be a non-empty character string
		obj <- name(object)
		if( !length(obj) || (length(obj)>1L || obj=='') )
			return(str_c("Slot 'name' must be a single non-empty character string [", obj, ']'))
		TRUE
	}
)



setGeneric('name', function(object, ...) standardGeneric('name'))

setMethod('name', signature(object='Strategy'),
	function(object, all=FALSE){
		n <- slot(object, 'name')
		if( length(n) && !all ) n[1L] else n
	}
)

setGeneric('name<-', function(object, ..., value) standardGeneric('name<-'))

setReplaceMethod('name', signature(object='Strategy', value='character'),
	function(object, value){
		slot(object, 'name') <- value
		validObject(object)
		object
	}
)

defaultArgument <- function(name, object, value, force=FALSE){

	# taken from methods::hasArg
	aname <- as.character(substitute(name))
	miss <- eval(substitute(missing(name)), sys.frame(sys.parent()))
	defaults <- attr(object, 'defaults')

	if( !miss && !force ) eval(substitute(name), sys.frame(sys.parent()))
	else if( aname %in% names(defaults) ) defaults[[aname]] 
	else value
}


setClass('KINOMOStrategy'
	, representation(
				objective = '.functionSlot' # the objective function used to compute the error (defined by name or function)
				, model = 'character' # KINOMO model to use
				, mixed = 'logical' # can the input data be negative?
	)
	, prototype=prototype(objective='euclidean', model='KINOMOstd', mixed=FALSE)
	, validity=function(object){
		
		# slot 'objective' must either be a non-empty character string or a function
		obj <- objective(object)
		if( is.character(obj) && obj == '' )
			return("Slot 'objective' must either be a non-empty character string or a function definition.")
			
		# slot 'model' must be the name of a class that extends class 'KINOMO'
		obj <- modelname(object)
		if( !is.character(obj) )
			return("Slot 'model' must be a character vector")
		if( any(inv <- !sapply(obj, isKINOMOclass)) )
			return(paste("Slot 'model' must contain only names of a class that extends class 'KINOMO' [failure on class(es) "
					, paste( paste("'", obj[inv], "'", sep=''), collapse=', ')  
					,"]"
					, sep=''))
		
		# slot 'mixed' must be a single logical		
		obj <- slot(object, 'mixed')
		if( length(obj) != 1 )
			return( paste("Slot 'mixed' must be a single logical [length=", length(obj), "]", sep='') )
	}
	, contains = c('VIRTUAL', 'Strategy')
)

#' @export
#' @rdname KINOMOStrategy-class
setMethod('show', 'KINOMOStrategy',
		function(object){			
			cat('<object of class: ', class(object), ">\n", sep='')
			cat(" name: ", name(object), " [", packageSlot(object), "]\n", sep='')
			svalue <- objective(object)
			svalue <- if( is.function(svalue) ) str_args(svalue, exdent=10) else paste("'", svalue,"'", sep='')
			cat(" objective:", svalue, "\n")
			cat(" model:", modelname(object), "\n")
			if( length(object@defaults) > 0L ){
				cat(" defaults:", str_desc(object@defaults, exdent=10L), "\n")
			}
			return(invisible())
		}
)

# Coerce method for 'KINOMOStrategy' objects into 'character': give the main name
setAs('KINOMOStrategy', 'character'
	, def = function(from) name(from)	
) 


setGeneric('KINOMOStrategy', function(name, method, ...) standardGeneric('KINOMOStrategy') )

setMethod('KINOMOStrategy', signature(name='character', method='function'), 
		function(name, method, ...){
			
			# build a KINOMOStrategyFunction object on the fly to wrap function 'method'
			KINOMOStrategy(name=name, algorithm=method, ...)
			
		}
)

#' Creates an \code{KINOMOStrategy} object based on a template object (Constructor-Copy).
setMethod('KINOMOStrategy', signature(name='character', method='KINOMOStrategy'), 
		function(name, method, ...){
			
			package <- topns_name()
			# build an KINOMOStrategy object based on template object
			strategy <- new(class(method), method, name=name, ..., package=package)
			
			# valid the new strategy
			validObject(strategy)
			
			# add trace of inheritance from parent KINOMO algorithm
			attr(strategy, 'parent') <- name(method)[1]
			
			# return new object
			strategy
		}
)

#' Creates an \code{KINOMOStrategy} based on a template object (Constructor-Copy), 
#' in particular it uses the \strong{same} name.
setMethod('KINOMOStrategy', signature(name='KINOMOStrategy', method='missing'), 
		function(name, method, ...){
			
			# do not change the object if single argument
			if( nargs() == 1L ) return(name)
			
			# use the name as a key
			# NB: need special trick to avoid conflict between argument and function 
			mname <- match.fun('name')(name)

			KINOMOStrategy(name=mname, method=name, ...)
		}
)


setMethod('KINOMOStrategy', signature(name='missing', method='character'), 
		function(name, method, ...){
			KINOMOStrategy(KINOMOAlgorithm(method, exact=TRUE), ...)
		}
)



setMethod('KINOMOStrategy', signature(name='NULL', method='KINOMOStrategy'), 
		function(name, method, ...){
			
			# use the name as a key
			# NB: need special trick to avoid conflict between argument and function 
			mname <- match.fun('name')(method)
			mname <- basename(tempfile(str_c(mname, '_')))
			
			KINOMOStrategy(name=mname, method=method, ...)
		}
)


setMethod('KINOMOStrategy', signature(name='character', method='character'), 
		function(name, method, ...){
			KINOMOStrategy(name=name, method=KINOMOAlgorithm(method, exact=TRUE), ...) 
		}
)

setMethod('KINOMOStrategy', signature(name='NULL', method='character'), 
		function(name, method, ...){
			KINOMOStrategy(NULL, method=KINOMOAlgorithm(method, exact=TRUE), ...) 
		}
)

setMethod('KINOMOStrategy', signature(name='character', method='missing'), 
		function(name, method, ...){
			
			package <- topns_name()
			# check iterative strategy
			if( hasArg2('Update') ){ # create a new KINOMOStrategyIterative object
				new('KINOMOStrategyIterative', name=name, ..., package=package)
			}else if( hasArg2('algorithm') ){
				new('KINOMOStrategyFunction', name=name, ..., package=package)
			}else{
				stop('KINOMOStrategy - Could not infer the type of KINOMO strategy to instantiate.')
			}
			
		}
)


setMethod('run', signature(object='KINOMOStrategy', y='matrix', x='KINOMOfit'),
	function(object, y, x, ...){
		stop("KINOMOStrategy::run is a pure virtual method that should be overloaded in class '", class(object),"'.")
	}
)

setMethod('run', signature(object='KINOMOStrategy', y='matrix', x='KINOMO'),
	function(object, y, x, ...){
		run(object, y, KINOMOfit(fit=x, seed='none', method=name(object)), ...)
	}
)


setMethod('deviance', 'KINOMOStrategy',
	function(object, x, y, ...){
		
		obj.fun <- slot(object, 'objective')
		
		# return the distance computed using the strategy's objective function
		if( !is.function(obj.fun) )
			deviance(x, y, method=obj.fun, ...)
		else # directly compute the objective function
			obj.fun(x, y, ...)
		
	}
)
		

setMethod('objective', 'KINOMOStrategy',
	function(object){
		slot(object, 'objective')
	}
)

setReplaceMethod('objective', signature(object='KINOMOStrategy', value='character'),
	function(object, value){
		#TODO: test for the existence of objective method
		slot(object, 'objective') <- value
		validObject(object)
		object
	}
)

setReplaceMethod('objective', signature(object='KINOMOStrategy', value='function'),
	function(object, value){
		slot(object, 'objective') <- value
		validObject(object)
		object
	}
)


setMethod('modelname', signature(object='KINOMOStrategy'),
	function(object){
		slot(object, 'model')
	}
)

is.mixed <-	function(object){
	return( slot(object, 'mixed') )
}


KINOMOFormals <- function(x, ...){
	UseMethod('KINOMOFormals')
}

#' @export
KINOMOFormals.character <- function(x, ...){
	s <- KINOMOAlgorithm(x)
	KINOMOFormals(s, ...)
}

#' @export
KINOMOFormals.KINOMOStrategy <- function(x, ...){
	m <- getMethod('run', signature(object='KINOMOStrategy', y='matrix', x='KINOMOfit'))
	args <- allFormals(m)
	# prepend registered default arguments
	expand_list(x@defaults, args)
}


KINOMOArgs <- function(x){
	args(KINOMOWrapper(x))
}
