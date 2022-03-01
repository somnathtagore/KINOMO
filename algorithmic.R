

#' @name algorithmic-KINOMO
#' @rdname algorithmic
NULL


#' @rdname algorithmic
setGeneric('algorithm<-', function(object, ..., value) standardGeneric('algorithm<-') )


#' @rdname algorithmic
setGeneric('seeding', function(object, ...) standardGeneric('seeding') )

setGeneric('seeding<-', function(object, ..., value) standardGeneric('seeding<-') )


#' @rdname algorithmic
setGeneric('niter', function(object, ...) standardGeneric('niter'))

#' @export
setGeneric('niter<-', function(object, ..., value) standardGeneric('niter<-'))



#' @rdname algorithmic
setGeneric('nrun', function(object, ...) standardGeneric('nrun') )

setMethod('nrun', 'ANY', 
	function(object){
		attr(object, 'nrun')
	}
)


#' @rdname algorithmic 
setGeneric('objective', function(object, ...) standardGeneric('objective'))

#' @export
#' @rdname algorithmic
setGeneric('objective<-', function(object, ..., value) standardGeneric('objective<-'))


setGeneric('runtime', function(object, ...) standardGeneric('runtime') )

#' @rdname algorithmic
setGeneric('runtime.all', function(object, ...) standardGeneric('runtime.all') )


#' @rdname algorithmic
setGeneric('seqtime', function(object, ...) standardGeneric('seqtime') )


#' @export
setGeneric('modelname', function(object, ...) standardGeneric('modelname'))


setMethod('modelname', 'ANY', 
	function(object)
	{
		as.character(class(object))
	}
)


#' @rdname algorithmic
setGeneric('run', function(object, y, x, ...) standardGeneric('run'))


setGeneric('logs', function(object, ...) standardGeneric('logs'))

#' It returns \code{NULL} if no logging data was found.
setMethod('logs', 'ANY', 
	function(object)
	{
		res <- attr(object, 'logs')
		if( !is.null(res) ) res 
		else if( is.list(object) ) object$logs
	}
)


#' @rdname algorithmic
setGeneric('compare', function(object, ...) standardGeneric('compare') )
