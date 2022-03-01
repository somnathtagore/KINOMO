#' @include registry.R
#' @include KINOMOStrategy-class.R
NULL


#'
setClass('KINOMOSeed'
		, representation(
			method = 'function' # the method actual definition
		)
		, contains = 'Strategy'
)

#' Show method for objects of class \code{KINOMOSeed} 
setMethod('show', 'KINOMOSeed',
		function(object){			
			cat('<object of class: ', class(object), ">\n")
			cat("name:\t", name(object), "\n")
			svalue <- algorithm(object)
			svalue <- if( is.function(svalue) ) '<function>' else paste("'", svalue,"'", sep='')
			cat("method:\t", svalue, "\n")
			return(invisible())
		}
)

#' Returns the workhorse function of the seeding method described by \code{object}. 
setMethod('algorithm', signature(object='KINOMOSeed'),
	function(object){						
		slot(object, 'method')
	}
)
#' Sets the workhorse function of the seeding method described by \code{object}.
setReplaceMethod('algorithm', signature(object='KINOMOSeed', value='function'),
	function(object, value){
		slot(object, 'method') <- value
		validObject(object)
		object
	}
)


setGeneric('KINOMOSeed', function(key, method, ...) standardGeneric('KINOMOSeed') )
#' Default method simply calls \code{\link{new}} with the same arguments. 
setMethod('KINOMOSeed', signature(key='character', method='ANY'), 
		function(key, method, ...){
			# wrap function method into a new KINOMOSeed object
			new('KINOMOSeed', name=key, method=method, ..., package=topns_name())
		}
)

#' Creates an \code{KINOMOSeed} based on a template object (Constructor-Copy), 
#' in particular it uses the \strong{same} name.
setMethod('KINOMOSeed', signature(key='KINOMOSeed', method='ANY'), 
		function(key, method, ...){
			
			# do not change the object if single argument
			if( nargs() == 1L ) return(key)
			
			# build an object based on template object
			new(class(method), key, method=method, ..., package=topns_name())
			
		}
)

