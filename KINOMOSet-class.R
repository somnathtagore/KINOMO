#' @include KINOMOfit-class.R
#' @include heatmaps.R
NULL



#' isKINOMOfit(list(res, resm, 'not a result'), recursive=FALSE)
#' 	 
isKINOMOfit <- function(object, recursive=TRUE){
	res <- is(object, 'KINOMOfit') || is(object, 'KINOMOfitX')
	# if the object is not a KINOMO result: apply to each element if a list (only in recursive mode)
	if( !res  && recursive && is.list(object) )
		sapply(object, isKINOMOfit)
	else
		res
}


#'   
setClass('KINOMOList'
		, representation(
			runtime='proc_time'
		)
		, contains='namedList'
		, validity=function(object){
			
			# the list must only contains KINOMOfit objects of the same dimensions
			ok <- isKINOMOfit(object)
			if( !is.logical(ok) )
				return("Could not validate elements in list: input is probably a complex structure of lists.")
			pb <- which(!ok)
			if( length(pb) ){
				return(paste("invalid class for element(s)"
							, str_out(i)
							, "of input list [all elements must be fitted KINOMO models]"))	
			}			
		}
)

#' Show method for objects of class \code{KINOMOList}
#' @export
setMethod('show', 'KINOMOList', 
	function(object)
	{
		cat("<Object of class:", class(object), ">\n")
		cat("Length:", length(object), "\n")
		if( length(object) > 0 ) cat("Method(s):", algorithm(object, string=TRUE), "\n")
		# show totaltime if present
		tt <- runtime(object)
		if( length(tt) > 0 ){
			cat("Total timing:\n"); show(tt);
		}
	}
)


#' This argument is forced to \code{TRUE} when \code{string=TRUE}.
#' 
setMethod('algorithm', 'KINOMOList', 
	function(object, string=FALSE, unique=TRUE){
		l <- length(object)
		
		if( string ) unique <- TRUE
		
		if( l == 0 ) NULL
		else if( l == 1 ) algorithm(object[[1]])
		else{	
			# build the vector of the algorithm names (with no repeat)
			m <- sapply(object, algorithm)
			if( unique ) m <- unique(m)
			if( string ) m <- paste(m, collapse=', ')
			m
		}
	} 
)

.seqtime <- function(object){
	
	if( length(object) == 0 ) return(NULL)
	# sum up the time across the runs
	t.mat <- sapply(object, function(x){
		if( is(x, 'KINOMOfitXn') ) runtime.all(x)
		else runtime(x)
	})
	res <- rowSums(t.mat)
	class(res) <- 'proc_time'
	res
}


#' It returns \code{NULL} on an empty object.
setMethod('seqtime', 'KINOMOList', 
	function(object){	
		if( length(object) == 0 ) return(NULL)
		# sum up the time across the runs
		.seqtime(object)
	}
)


#' the fits in \code{object}.
setMethod('runtime', 'KINOMOList', 
	function(object, all=FALSE){
		if( !all ){
			t <- slot(object, 'runtime')
			if( length(t)==0 ) seqtime(object) else t
		}else
			sapply(object, runtime)
	}
)

as.KINOMOList <- function(..., unlist=FALSE){
	arg.l <- list(...)
	if( length(arg.l) == 1L && is.list(arg.l[[1]]) && !is(arg.l[[1]], 'KINOMOfitX') )
		arg.l <- arg.l[[1]]
	
	# unlist if required
	if( unlist )
		arg.l <- unlist(arg.l)
	
	# create a KINOMOList object from the input list 
	new('KINOMOList', arg.l)
}
			


#' res
#' 
#' # plot a heatmap of the consensus matrix
#' \dontrun{ consensusmap(res) }
#'  
setClass('KINOMOfitX'
		, representation(
				runtime.all = 'proc_time' # running time to perform all the KINOMO runs
		)
		, contains='VIRTUAL'		
)

#' Returns the CPU time required to compute all the KINOMO runs.
#' It returns \code{NULL} if no CPU data is available. 
setMethod('runtime.all', 'KINOMOfitX', 
		function(object){	
			t <- slot(object, 'runtime.all')
			if( length(t) > 0 ) t else NULL
		}
)

#
' greater than one, while only the result of the best run is stored in 
#' the object (cf. option \code{'k'} in method \code{\link{KINOMO}}).
setMethod('nrun', 'KINOMOfitX', 
	function(object){	
		stop("KINOMO::KINOMOfitX - missing definition for pure virtual method 'nrun' in class '", class(object), "'")
	}
)

#' This method always returns 1, since an \code{KINOMOfit} object is obtained 
#' from a single KINOMO run.  
setMethod('nrun', 'KINOMOfit', 
	function(object){
		1L
	}
)


#' @rdname connectivity
#' @export
setGeneric('consensus', function(object, ...) standardGeneric('consensus') )
#' Pure virtual method defined to ensure \code{consensus} is defined for sub-classes of \code{KINOMOfitX}.
#' It throws an error if called.
setMethod('consensus', 'KINOMOfitX', 
	function(object, ...){	
		stop("KINOMO::KINOMOfitX - missing definition for pure virtual method 'consensus' in class '", class(object), "'")
	}
)

#' This method is provided for completeness and is identical to 
#' \code{\link{connectivity}}, and returns the connectivity matrix, 
#' which, in the case of a single KINOMO model, is also the consensus matrix.
setMethod('consensus', 'KINOMO', 
	function(object, ...){
		connectivity(object, ...)
	}
)


#' @inline
#' @export
setGeneric('consensushc', function(object, ...) standardGeneric('consensushc'))


#' Default value is \code{TRUE}.
setMethod('consensushc', 'matrix', 
	function(object, method='average', dendrogram=TRUE){
		
		# hierachical clustering based on the connectivity matrix
		hc <- hclust(as.dist(1-object), method=method)
		
		# convert into a dendrogram if requested
		if( dendrogram ) as.dendrogram(hc)
        else hc
	}
)
#' Compute the hierarchical clustering on the connectivity matrix of \code{object}.
setMethod('consensushc', 'KINOMO', 
	function(object, ...){		
		# hierachical clustering based on the connectivity matrix
		consensushc(connectivity(object), ...)		
	}
)


#' 
setMethod('consensushc', 'KINOMOfitX', 
	function(object, what=c('consensus', 'fit'), ...){
		
		what <- match.arg(what)
		if( what == 'consensus' ){
			# hierachical clustering on the consensus matrix
			consensushc(consensus(object), ...)
						
		}else if( what == 'fit' )
			consensushc(fit(object), ...)
		
	}
)


setMethod('predict', signature(object='KINOMOfitX'),
	function(object, what=c('columns', 'rows', 'samples', 'features', 'consensus', 'chc'), dmatrix = FALSE, ...){
		# determine which prediction to do
		what <- match.arg(what)
		res <- if( what %in% c('consensus', 'chc') ){
			# build the tree from consensus matrix
			h <- consensushc(object, what='consensus', dendrogram=FALSE)
			# extract membership from the tree
			cl <- cutree(h, k=nbasis(object))
			
			# rename the cluster ids in the case of a consensus map
			if( what != 'chc' ){
                dr <- as.dendrogram(h)
                o <- order.dendrogram(reorder(dr, rowMeans(consensus(object), na.rm=TRUE)))
				cl <- setNames(match(cl, unique(cl[o])), names(cl))
            }
            
			res <- as.factor(cl)
            # add dissimilarity matrix if requested
            if( dmatrix ){
                attr(res, 'dmatrix') <- 1 - consensus(object) 
            }
            if( what != 'chc' ) attr(res, 'iOrd') <- o
            
            # return
            res
		}
		else predict(fit(object), what=what, ..., dmatrix = dmatrix)
        attr(res, 'what') <- what
        res
	}
)

#' Returns the model object that achieves the lowest residual approximation 
#' error across all the runs.
#' 
#' It is a pure virtual method defined to ensure \code{fit} is defined 
#' for sub-classes of \code{KINOMOfitX}, which throws an error if called.
setMethod('fit', 'KINOMOfitX', 
	function(object){	
		stop("KINOMO::KINOMOfitX - missing definition for pure virtual method 'fit' in class '", class(object), "'")
	}
)



#' for sub-classes of \code{KINOMOfitX}, which throws an error if called.
setMethod('minfit', 'KINOMOfitX',
	function(object){
		stop("KINOMO::KINOMOfitX - missing definition for pure virtual method 'minfit' in class '", class(object), "'")
	}
)

#' Show method for objects of class \code{KINOMOfitX}
#' @export
setMethod('show', 'KINOMOfitX', 
		function(object){
			cat("<Object of class:", class(object), ">\n")
			# name of the algorithm
			cat("  Method:", algorithm(object), "\n")
			# number of runs
			cat("  Runs: ", nrun(object),"\n");
			# initial state
			cat("  RNG:\n  ", RNGstr(getRNG1(object)),"\n");
			if( nrun(object) > 0 ){
				# show total timing			
				cat("  Total timing:\n"); show(runtime.all(object));
			}
		}
)


setGeneric('getRNG1', package='rngtools')


setMethod('getRNG1', signature(object='KINOMOfitX'),
	function(object){
		stop("KINOMO::getRNG1(", class(object), ") - Unimplemented pure virtual method: could not extract initial RNG settings.")
	}
)

#' Compares two KINOMO models when at least one comes from multiple KINOMO runs. 
setMethod('KINOMO.equal', signature(x='KINOMOfitX', y='KINOMO'), 
	function(x, y, ...){
		KINOMO.equal(fit(x), y, ...)
	}
)
#' Compares two KINOMO models when at least one comes from multiple KINOMO runs.
setMethod('KINOMO.equal', signature(x='KINOMO', y='KINOMOfitX'), 
		function(x, y, ...){
			KINOMO.equal(x, fit(y), ...)
		}
)

#' Returns the residuals achieved by the best fit object, i.e. the lowest 
#' residual approximation error achieved across all KINOMO runs.
setMethod('residuals', signature(object='KINOMOfitX'), 
		function(object, ...){
			residuals(minfit(object), ...)
		}
)
#' Returns the deviance achieved by the best fit object, i.e. the lowest 
#' deviance achieved across all KINOMO runs.
setMethod('deviance', signature(object='KINOMOfitX'), 
		function(object, ...){
			deviance(minfit(object), ...)
		}
)




#' 
setClass('KINOMOfitX1'
	, representation(
			#fit = 'KINOMOfit' # holds the best fit from all the runs
			consensus = 'matrix' # average connectivity matrix of all the KINOMO runs
			, nrun = 'integer'
			, rng1 = 'ANY'
	)
	, contains=c('KINOMOfitX', 'KINOMOfit')
	, prototype=prototype(
			consensus =	matrix(as.numeric(NA),0,0)
			, nrun = as.integer(0)
	)
)



#' Show method for objects of class \code{KINOMOfitX1}
#' @export
setMethod('show', 'KINOMOfitX1', 
	function(object){
		callNextMethod(object)
		
		# show details of the best fit
		#cat(" # Best fit:\n  ")
		#s <- capture.output(show(fit(object)))
		#cat(s, sep="\n  |")
	}
)

#' Returns the number of KINOMO runs performed, amongst which \code{object} was 
#' selected as the best fit.
setMethod('nrun', 'KINOMOfitX1', 
	function(object){
		slot(object,'nrun')
	}
)


#' The result is the matrix stored in slot \sQuote{consensus}.
#' This method returns \code{NULL} if the consensus matrix is empty.
setMethod('consensus', signature(object='KINOMOfitX1'), 
	function(object, no.attrib = FALSE){
		
		C <- slot(object, 'consensus')
		if( length(C) > 0 ){
			if( !no.attrib ){
				class(C) <- c(class(C), 'KINOMO.consensus')
				attr(C, 'nrun') <- nrun(object)
				attr(C, 'nbasis') <- nbasis(object)
			}
			C
		}else NULL
		
	}
)


#' Since \code{KINOMOfitX1} objects only hold the best fit, this method simply 
#' returns \code{object} coerced into an \code{KINOMOfit} object.
setMethod('minfit', 'KINOMOfitX1',
	function(object){	
		# coerce the object into a KINOMOfit object
		as(object, 'KINOMOfit')
	}
)


#' \sQuote{fit}.
setMethod('fit', signature(object='KINOMOfitX1'),
	function(object){
		slot(object, 'fit')
	}
)

#' Returns the RNG settings used to compute the first of all KINOMO runs, amongst
#' which \code{object} was selected as the best fit.
setMethod('getRNG1', signature(object='KINOMOfitX1'),
	function(object){
		object@rng1
	}
)

#' Compares the KINOMO models fitted by multiple runs, that only kept the best fits.
setMethod('KINOMO.equal', signature(x='KINOMOfitX1', y='KINOMOfitX1'), 
	function(x, y, ...){
		KINOMO.equal(fit(x), fit(y), ...)
	}
)




#' 
setClass('KINOMOfitXn'
	, contains=c('KINOMOfitX', 'list')
	, validity=function(object){
				
		# the list must only contains KINOMOfit objects of the same dimensions
		ref.dim <- NULL
		ref.algo <- NULL
		for(i in seq_along(object)){
			# check class of the element
			item <- object[[i]]
			if( !(is(item, 'KINOMOfit') && !is(item, 'KINOMOfitX')) )
				return(paste("invalid class for element", i, "of input list [all elements must be a KINOMOfit object]"))
						
			# check dimensions
			if( is.null(ref.dim) ) ref.dim <- dim(item)	
			if( !identical(ref.dim, dim(item)) )
				return(paste("invalid dimension for element", i, "of input list [all elements must have the same dimensions]"))
			
			# check algorithm names
			if( is.null(ref.algo) ) ref.algo <- algorithm(item)	
			if( !identical(ref.algo, algorithm(item)) )
				return(paste("invalid algorithm for element", i, "of input list [all elements must result from the same algorithm]"))
			
		}
		
	}
)


 
#' Show method for objects of class \code{KINOMOfitXn}
#' @export
setMethod('show', 'KINOMOfitXn', 
	function(object){
		callNextMethod(object)
		
		# if the object is not empty and slot runtime.all is not null then show
		# the sequential time, as it might be different from runtime.all
		if( length(object) > 0 && !is.null(runtime.all(object, null=TRUE)) ){
			# show total sequential timing
			cat("  Sequential timing:\n"); show(seqtime(object));
		}
	}
)


#' This method returns \code{NULL} if the object is empty.  
setMethod('nbasis', signature(x='KINOMOfitXn'), 
	function(x, ...){
		if( length(x) == 0 ) return(NULL)
		return( nbasis(x[[1]]) )
	}
)


#' @rdname dims
setMethod('dim', signature(x='KINOMOfitXn'), 
	function(x){
		if( length(x) == 0 ) return(NULL)
		return( dim(x[[1L]]) )
	}
)

#' Returns the coefficient matrix of the best fit amongst all the fits stored in 
#' \code{object}.
#' It is a shortcut for \code{coef(fit(object))}.  
setMethod('coef', signature(object='KINOMOfitXn'), 
	function(object, ...){
		coef(fit(object), ...)
	}
)

#' Returns the basis matrix of the best fit amongst all the fits stored in 
#' \code{object}.
#' It is a shortcut for \code{basis(fit(object))}.
setMethod('basis', signature(object='KINOMOfitXn'), 
	function(object, ...){
		basis(fit(object), ...)
	}
)

#' Method for multiple KINOMO fit objects, which returns the indexes of fixed basis 
#' terms from the best fitted model.
setMethod('ibterms', 'KINOMOfitX', 
	function(object){
		ibterms(fit(object))
	}
)
#' Method for multiple KINOMO fit objects, which returns the indexes of fixed 
#' coefficient terms from the best fitted model.
setMethod('icterms', 'KINOMOfit', 
	function(object){
		icterms(fit(object))
	}
)


#' Returns the number of runs performed to compute the fits stored in the list 
#' (i.e. the length of the list itself).
setMethod('nrun', 'KINOMOfitXn', 
	function(object){
		length(object)
	}
)

#' Returns the name of the common KINOMO algorithm used to compute all fits 
#' stored in \code{object}
#' 
#' Since all fits are computed with the same algorithm, this method returns the 
#' name of algorithm that computed the first fit.
#' It returns \code{NULL} if the object is empty.
setMethod('algorithm', 'KINOMOfitXn', 
	function(object){
		if( length(object) == 0 ) return(NULL)
		return( algorithm(object[[1]]) )
	} 
)
#' Returns the name of the common seeding method used the computation of all fits 
#' stored in \code{object}
#' 
#' Since all fits are seeded using the same method, this method returns the 
#' name of the seeding method used for the first fit.
#' It returns \code{NULL} if the object is empty.
setMethod('seeding', 'KINOMOfitXn',
	function(object){
		if( length(object) == 0 ) return(NULL)
		return( seeding(object[[1]]) )
	}
)
#' Returns the common type KINOMO model of all fits stored in \code{object}
#' 
#' Since all fits are from the same KINOMO model, this method returns the 
#' model type of the first fit.
#' It returns \code{NULL} if the object is empty.
setMethod('modelname', signature(object='KINOMOfitXn'), 
	function(object){
		if( length(object) == 0 ) return(NULL)
		return( modelname(object[[1]]) )
	}
)

#' Returns the CPU time that would be required to sequentially compute all KINOMO 
#' fits stored in \code{object}.
#' 
#' This method calls the function \code{runtime} on each fit and sum up the 
#' results.
#' It returns \code{NULL} on an empty object.
setMethod('seqtime', 'KINOMOfitXn', 
		function(object){	
			if( length(object) == 0 ) return(NULL)
			# sum up the time across the runs
			.seqtime(object)
		}
)



#' 
setMethod('runtime.all', 'KINOMOfitXn', 
	function(object, null=FALSE, warning=TRUE){
		
		if( length(object) == 0 ) return(NULL)
		stored.time <- slot(object, 'runtime.all')
		# if there is some time stored, return it
		if( length(stored.time) > 0 ) stored.time
		else if( null ) NULL
		else{
			if( warning )
				warning("KINOMOfitXn::runtime.all - computation time data not available [sequential time was used instead]")
			seqtime(object) # otherwise total sequential time
		}
		
	}
)


#'  
setMethod('minfit', 'KINOMOfitXn',
	function(object){
		
		b <- which.best(object, deviance)
		# test for length 0
		if( length(b) == 0 ) return(NULL)
				
		# return the run with the lower
		object[[ b ]]
	}
)



#' @rdname advanced
which.best <- function(object, FUN=deviance, ...){
	
	# test for length 0
	if( length(object) == 0 ) 
		return(integer())
	
	# retrieve the measure for each run
	e <- sapply(object, FUN, ...)
	
	# return the run with the lower
	which.min(e)
}


#' This method throws an error if the object is empty.
setMethod('getRNG1', signature(object='KINOMOfitXn'),
	function(object){
		if( length(object) == 0 )
			stop("KINOMO::getRNG1 - Could not extract RNG data from empty object [class:", class(object), "]")
		
		getRNG(object[[1]])
	}
)

#' @inline
#' @rdname RNG
#' @export
setGeneric('.getRNG', package='rngtools')

#' Returns the RNG settings used for the best fit.
#' 
#' This method throws an error if the object is empty.
setMethod('.getRNG', signature(object='KINOMOfitXn'),
	function(object, ...){
		if( length(object) == 0 )
			stop("KINOMO::getRNG - Could not extract RNG data from empty object [class:", class(object), "]")
		
		getRNG(minfit(object), ...)
	}
)


#' Returns the best KINOMO fit object amongst all the fits stored in \code{object},
#' i.e. the fit that achieves the lowest estimation residuals.
setMethod('fit', signature(object='KINOMOfitXn'),
	function(object){
		fit( minfit(object) )
	}
)


#' @inline
setMethod('KINOMO.equal', signature(x='list', y='list'), 
	function(x, y, ..., all=FALSE, vector=FALSE){
		if( !all )
			KINOMO.equal(x[[ which.best(x) ]], y[[ which.best(y) ]], ...)
		else{
			if( length(x) != length(y) )
				FALSE
			else
				res <- mapply(function(a,b,...) isTRUE(KINOMO.equal(a,b,...)), x, y, MoreArgs=list(...))
				if( !vector )
					res <- all( res )
				res
		}
	}
)

#' Compare all elements in \code{x} to \code{x[[1]]}.
setMethod('KINOMO.equal', signature(x='list', y='missing'), 
	function(x, y, ...){
		
		if( length(x) == 0L ){
			warning("Empty list argument `x`: returning NA")
			return(NA)
		}
		if( length(x) == 1L ){
			warning("Only one element in list argument `x`: returning TRUE")
			return(TRUE)
		}
		for( a in x ){
			if( !KINOMO.equal(x[[1]], a, ...) ) return(FALSE)
		}
		return(TRUE)
	}
)


#' @aliases plot.KINOMO.consensus
setMethod('consensus', signature(object='KINOMOfitXn'), 
	function(object, ..., no.attrib = FALSE){
		if( length(object) == 0 ) return(NULL)
		
		# init empty consensus matrix
		con <- matrix(0, ncol(object), ncol(object))
		# name the rows and columns appropriately: use the sample names of the first fit
		dimnames(con) <- list(colnames(object[[1]]), colnames(object[[1]]))
		
		# compute mean connectivity matrix
		sapply(object 
				, function(x, ...){
					con <<- con + connectivity(x, ..., no.attrib = TRUE)
					NULL
				}
				, ...
		)
		con <- con / nrun(object)
				
		# return result
		if( !no.attrib ){
			class(con) <- c(class(con), 'KINOMO.consensus')
			attr(con, 'nrun') <- nrun(object)
			attr(con, 'nbasis') <- nbasis(object)
		}
		con
	}
)

#' @method plot KINOMO.consensus
#' @export 
plot.KINOMO.consensus <- function(x, ...){
	consensusmap(x, ...)
}


#' @export 
setGeneric('dispersion', function(object, ...) standardGeneric('dispersion') )
#' Workhorse method that computes the dispersion on a given matrix.
setMethod('dispersion', 'matrix', 
	function(object, ...){
		stopifnot( nrow(object) == ncol(object) )
		sum( 4 * (object-1/2)^2 ) / nrow(object)^2
	}
)
#' Computes the dispersion on the consensus matrix obtained from multiple KINOMO
#' runs. 
setMethod('dispersion', 'KINOMOfitX', 
	function(object, ...){
		dispersion(consensus(object), ...)
	}
)


#' @keywords internal
setGeneric('KINOMOfitX', function(object, ...) standardGeneric('KINOMOfitX') )
#' Create an \code{KINOMOfitX} object from a list of fits.


#' 
setMethod('KINOMOfitX', 'list',
	function(object, ..., .merge=FALSE){
		
		if( length(object) == 0 )
			return(new('KINOMOfitXn'))
		else if( is(object, 'KINOMOfitXn') && !.merge)
			return(object)
		
		# retrieve the extra arguments
		extra <- list(...)
				
		# if runtime.all is provided: be sure it's of the right class
		tt <- extra$runtime.all
		compute.tt <- TRUE
		if( !is.null(tt) ){
			if( !is(tt, 'proc_time') ){
				if( !is.numeric(tt) || length(tt) != 5 )
					stop("KINOMO::KINOMOfitX - invalid value for 'runtime.all' [5-length numeric expected]")
				class(extra$runtime.all) <- 'proc_time'
			}
			compute.tt <- FALSE
		}else{
			extra$runtime.all <- rep(0,5)
			class(extra$runtime.all) <- 'proc_time'
		}
		
		# check validity and aggregate if required
		ref.algo <- NULL
		ref.class <- NULL
		nrun <- 0
		lapply( seq_along(object)
			, function(i){
				item <- object[[i]]
				
				# check the type of each element				
				if( !(is(item, 'KINOMOfitX') || is(item, 'KINOMOfit')) )
					stop("KINOMO::KINOMOfitX - invalid class for element ", i, " of input list [all elements must be KINOMOfit or KINOMOfitX objects]")
				
				# check that all elements result from the same algorithm
				if( is.null(ref.algo) ) ref.algo <<- algorithm(item)
				if( !identical(algorithm(item), ref.algo) )
					stop("KINOMO::KINOMOfitX - invalid algorithm for element ", i, " of input list [cannot join results from different algorithms]")
				
				# check if simple join is possible: only Ok if all elements are from the same class (KINOMOfit or KINOMOfitXn)
				if( length(ref.class) <= 1 ) ref.class <<- unique(c(ref.class, class(item)))
				
				# sum up the number of runs
				nrun <<- nrun + nrun(item)
				
				# compute total running time if necessary
				if( compute.tt )
					extra$runtime.all <<- extra$runtime.all + runtime.all(item)
				
			}
		)
		
		# force merging if the input list is hetergeneous or if it only contains KINOMOfitX1 objects
		if( length(ref.class) > 1 || ref.class == 'KINOMOfitX1' ){
			KINOMO.debug('KINOMOfitX', ".merge is forced to TRUE")
			.merge <- TRUE
		}
		
		# unpack all the KINOMOfit objects
		object.list <- unlist(object)
		KINOMO.debug('KINOMOfitX', "Number of fits to join = ", length(object.list))
					
		# one wants to keep only the best result
		if( .merge ){
			
			warning("KINOMO::KINOMOfitX - The method for merging lists is still in development")
			
			# set the total number of runs
			extra$nrun <- as.integer(nrun)			
									
			# consensus matrix
			if( !is.null(extra$consensus) )
				warning("KINOMO::KINOMOfitX - the value of 'consensus' was discarded as slot 'consensus' is computed internally")
			extra$consensus <- NULL
									
			consensus <- matrix(as.numeric(NA), 0, 0)
			best.res <- Inf		
			best.fit <- NULL
			sapply(object.list, function(x){
				if( !is(x, 'KINOMOfit') )
					stop("KINOMO::KINOMOfitX - all inner-elements of '",substitute(object),"' must inherit from class 'KINOMOfit'")
				
				# merge consensus matrices
				consensus <<- if( sum(dim(consensus)) == 0 ) nrun(x) * consensus(x)
							  else consensus + nrun(x) * consensus(x)
					  
				temp.res <- residuals(x)
				if( temp.res < best.res ){
					# keep best result
					best.fit <<- minfit(x)					
					best.res <<- temp.res
				}
			})
			# finalize consensus matrix
			consensus <- consensus/extra$nrun
			extra$consensus <- consensus 
									
			# return merged result
			return( do.call(KINOMOfitX, c(list(best.fit), extra)) )
		}
		else{
			# create a KINOMOfitXn object that holds the whole list			
			do.call('new', c(list('KINOMOfitXn', object.list), extra))
		}
	}
)
#' Creates an \code{KINOMOfitX1} object from a single fit.
#' This is used in \code{\link{KINOMO}} when only the best fit is kept in memory or 
#' on disk.
#'  
setMethod('KINOMOfitX', 'KINOMOfit',
		function(object, ...){
					
			extra <- list(...)
						
			# default value for nrun is 1 
			if( is.null(extra$nrun) ) extra$nrun = as.integer(1)
			
			# a consensus matrix is required (unless nrun is 1)
			if( is.null(extra$consensus) ){
				if( extra$nrun == 1 ) 
					extra$consensus <- connectivity(object)
				else
					stop("Slot 'consensus' is required to create a 'KINOMOfitX1' object where nrun > 1")				
			}
			
			# slot runtime.all is inferred if missing and nrun is 1
			if( is.null(extra$runtime.all) && extra$nrun == 1 )
				extra$runtime.all <- runtime(object)
			
			# create the KINOMOfitX1 object
			do.call('new', c(list('KINOMOfitX1', object), extra))
		}
)
#' Provides a way to aggregate \code{KINOMOfitXn} objects into an \code{KINOMOfitX1} 
#' object.
setMethod('KINOMOfitX', 'KINOMOfitX',
		function(object, ...){

			# nothing to do in the case of KINOMOfitX1 objects
			if( is(object, 'KINOMOfitX1') ) return(object)
			
			# retrieve extra arguments
			extra <- list(...)
			
			# take runtime.all from the object itself
			if( !is.null(extra$runtime.all) )
				warning("KINOMO::KINOMOfitX - argument 'runtime.all' was discarded as it is computed from argument 'object'")			
			extra$runtime.all <- runtime.all(object)						
			
			# create the KINOMOfitX1 object
			f <- selectMethod(KINOMOfitX, 'list')
			do.call(f, c(list(object), extra))
		}
)

#' Computes the best or mean purity across all KINOMO fits stored in \code{x}.
#' 
#' @param method a character string that specifies how the value is computed.
#' It may be either \code{'best'} or \code{'mean'} to compute the best or mean 
#' purity respectively.
#'  
#' @inline
setMethod('purity', signature(x='KINOMOfitXn', y='ANY'), 
	function(x, y, method='best', ...){
		c <- sapply(x, purity, y=y, ...)
		
		# aggregate the results if a method is provided
		if( is.null(method) ) c
		else aggregate.measure(c, method, decreasing=TRUE)		
	}
)

#' Computes the best or mean entropy across all KINOMO fits stored in \code{x}.
#' 
#' @inline
setMethod('entropy', signature(x='KINOMOfitXn', y='ANY'), 
	function(x, y, method='best', ...){
		c <- sapply(x, entropy, y=y, ...)		
		
		# aggregate the results if a method is provided
		if( is.null(method) ) c
		else aggregate.measure(c, method)
	}
)


aggregate.measure <- function(measure, method=c('best', 'mean'), decreasing=FALSE){
	# aggregate the results
	method <- match.arg(method)
	res <- switch(method
			, mean = mean(measure)
			, best = if( decreasing ) max(measure) else min(measure)
	)
	
	# set the name to 
	names(res) <- method
	
	# return result
	res
} 



#'   
setMethod('summary', signature(object='KINOMOfitX'),
	function(object, ...){
		
		# compute summary measures for the best fit
		best.fit <- minfit(object)
		s <- summary(best.fit, ...)
		# get totaltime
		t <- runtime.all(object)		
		
		# replace cpu.all and nrun in the result (as these are set by the summary method of class KINOMOfit)
		s[c('cpu.all', 'nrun')] <- c(as.numeric(t['user.self']+t['user.child']), nrun(object))
		
		# compute cophenetic correlation coeff and dispersion
		C <- consensus(object)
		s <- c(s, cophenetic=cophcor(C), dispersion=dispersion(C))
		
        # compute mean consensus silhouette width
		si <- silhouette(object, what = 'consensus')
		s <- c(s, silhouette.consensus = if( !is_NA(si) ) summary(si)$avg.width else NA)
        
		# return result
		s
	}
)



#' @name compare-KINOMO
#' @rdname KINOMO-compare
NULL

.compare_KINOMO <- function(...){
	args <- list(...)
	
	iargs <-
	if( is.null(names(args)) ){
		names(args) <- rep("", length(args)) 
		seq(args)
	}else{
		iargs <- which(names(args)=='')
		if( length(iargs) != length(args) )
			iargs <- iargs[ iargs < which(names(args)!='')[1L] ]
		iargs
	}

	lfit <- args[iargs]
	lfit <- unlist(lfit, recursive=FALSE)
	
	# wrap up into an KINOMOList object
	object <- as.KINOMOList(lfit)
	do.call('summary', c(list(object), args[-iargs]))
}


#' 
setMethod('compare', signature(object='KINOMOfit'),
	function(object, ...){
		.compare_KINOMO(object, ...)
	}
)


setMethod('compare', signature(object='KINOMOfitXn'),
	function(object, ...){
		do.call(.compare_KINOMO, c(unlist(object), list(...)))
	}
)


#' # compare elements of a list
#' compare(list(res, res2), target=x)
setMethod('compare', signature(object='list'),
	function(object, ...){
		do.call(.compare_KINOMO, c(list(object), list(...)))
	}
)


#' 
#' @rdname KINOMO-compare
setMethod('summary', signature(object='KINOMOList'),
	function(object, sort.by=NULL, select=NULL, ...){
		
		if( length(object) == 0L ) return()
		
		# define the sorting schema for each criteria (TRUE for decreasing, FALSE for increasing)
		sorting.schema <- list(method=FALSE, seed=FALSE, rng=FALSE, metric=FALSE
							, residuals=FALSE, cpu=FALSE, purity=TRUE, nrun=FALSE, cpu.all=FALSE
							, cophenetic=TRUE, dispersion=TRUE #KINOMOfitX only
							, entropy=FALSE, sparseness.basis=TRUE, sparseness.coef=TRUE, rank=FALSE, rss=FALSE
							, niter=FALSE, evar=TRUE
                            , silhouette.coef = TRUE, silhouette.basis = TRUE
                            , silhouette.consensus = TRUE)
				
		# for each result compute the summary measures
		measure.matrix <- sapply(object, summary, ...)		
		
		# the results from 'summary' might not have the same length => generate NA where necessary
		if( is.list(measure.matrix) ){
			name.all <- unique(unlist(sapply(measure.matrix, names)))
			measure.matrix <- sapply(seq_along(measure.matrix),
				function(i){
					m <- measure.matrix[[i]][name.all]
					names(m) <- name.all
					m
				}
			)
		}
		
		# transpose the results so that methods are in lines, measures are in columns
		measure.matrix <- t(measure.matrix)	
		
		# set up the resulting data.frame
		methods <- sapply(object, function(x, ...){
					x <- minfit(x)
					m <- algorithm(x)
					s <- seeding(x) 
					svalue <- objective(x)
					svalue <- if( is.function(svalue) ) '<function>' else svalue
					c(method=m, seed=s, rng=RNGdigest(x), metric=svalue)
				}
		)
		methods <- t(methods)	
		res <- as.data.frame(methods, stringsAsFactors=FALSE)	
		
		# add the measures to the result		
		res <- cbind(res, measure.matrix)
		res$rng <- as.numeric(factor(res$rng))
				
		# sort according to the user's preference
		# ASSERT FOR DEV: all columns measure must have a defined sorting schema 
		#if( !all( no.schema <- is.element(colnames(res), names(sorting.schema))) ) 
		#	warning("ASSERT: missing sorting schema for criteria(e): ", paste(paste("'", colnames(res)[!no.schema], "'", sep=''), collapse=', '))

		if( !is.null(sort.by) ){			
			sorting.criteria <- intersect(colnames(res), names(sorting.schema))
			sort.by.ind <- pmatch(sort.by, sorting.criteria)
			if( is.na(sort.by.ind) )
				stop("KINOMO::summary[KINOMOList] : argument 'sort.by' must be NULL or partially match one of "
					, paste( paste("'", names(sorting.schema), "'", sep=''), collapse=', ')
					, call.=FALSE)
			sort.by <- sorting.criteria[sort.by.ind]
			res <- res[order(res[[sort.by]], decreasing=sorting.schema[[sort.by]]) , ]
			
			# add an attribute to the result to show the sorting criteria that was used
			attr(res, 'sort.by') <- sort.by
		}
		
		# limit the output to the required measures
		if( !is.null(select) || !missing(select) ){
			select.full <- match.arg(select, colnames(res), several.ok=TRUE)
			if( length(select.full) <  length(select) )
				stop("KINOMO::summary[KINOMOList] - the elements of argument 'select' must partially match one of "
					, paste(paste("'", colnames(res),"'", sep=''), collapse=', ')
					, call.=FALSE)
			res <- subset(res, select=select.full)
		}
				
		# return result
		res
	}
)


#' 
#' @rdname KINOMO-compare
setMethod('plot', signature(x='KINOMOList', y='missing'), 
	function(x, y, skip=-1L, ...){
		
		# retrieve normalized residuals tracks
		max.iter <- 0
		tracks <- lapply( x, 
				function(res){
					res <- minfit(res)
					t <- residuals(res, track=TRUE)
					# skip some residuals(s) if requested
					if( skip == -1L && !is.null(names(t)) ) t <- t[names(t)!='0'] # remove initial residual
					else if( skip > 0 ) t <- t[-(1:skip)]
					#print(t)
					# update max iteration
					max.iter <<- max(max.iter, as.numeric(names(t)))
					# return normalized track
					t/t[1]
				}
		)
		minT <- min(sapply(tracks, min))
		maxT <- max(sapply(tracks, max))
		
		#print(tracks)
		# create an empty plot
		# set default graphical parameters (those can be overriden by the user)
		params <- .set.list.defaults(list(...)
				, xlab='Iterations', ylab='Normalised objective values'
				, main='KINOMO Residuals')
		
		# setup the plot
		do.call('plot', 
				c(list(0, xlim=c(0,max.iter+100), ylim=c(minT, maxT)), col='#00000000'
				, params)
				)
		
		# add legend
		cols <- seq_along(tracks)
		legend('topright', legend=names(tracks), fill=cols
				, title='Algorithm')
		
		# plot each tracks		
		lapply( seq_along(tracks),
				function(i){
					t <- tracks[[i]]
					points(names(t), t, col=cols[i], type='p', cex=0.5)
					points(names(t), t, col=cols[i], type='l', lwd=1.4)
				})
		
		
		# return invisible
		return(invisible())
	}		
)

#' Deprecated method subsituted by \code{\link{consensusmap}}.
setMethod('metaHeatmap', signature(object='KINOMOfitX'),
		function(object, ...){
			# send deprecated warning
			.Deprecated('metaHeatmap', 'KINOMO', "Direct use of the S4-Method 'metaHeatmap' for 'KINOMOfitX' objects is deprecated, use 'consensusmap' instead.")

			# call the new function 'consmap'
			return( consensusmap(object, ...) )
			
		}
)



#' @export
setGeneric('consensusmap', function(object, ...) standardGeneric('consensusmap') )
#' Plots a heatmap of the consensus matrix obtained when fitting an KINOMO model with multiple runs. 
setMethod('consensusmap', 'KINOMOfitX', 
	function(object, annRow=NA, annCol=NA
			, tracks=c('basis:', 'consensus:', 'silhouette:')
			, main = 'Consensus matrix', info = FALSE
			, ...){
			
		# add side information if requested
		info <- if( isTRUE(info) ){
					paste("KINOMO model: '", modelname(object)
					, "'\nAlgorithm: '", algorithm(object)
					, "'\nbasis: ", nbasis(object)
					,"\nnrun: ", nrun(object), sep='')
				}else if( isFALSE(info) ) NULL
				else info
		
		x <- consensus(object)
		
		# process annotation tracks
		ptracks <- process_tracks(x, tracks, annRow, annCol)
		annRow <- ptracks$row 
		annCol <- ptracks$col
		# set special annotation handler
		ahandlers <- list(
			basis = function() predict(object)
			, consensus = function() predict(object, what='consensus')
			, silhouette = function(){
				si <- silhouette(object, what='consensus', order = NA)
				if( is_NA(si) ) NA
				else si[, 'sil_width']
			}
		)
		specialAnnotation(1L, ahandlers)
		specialAnnotation(2L, ahandlers)
		#
		
		consensusmap(x, ..., annRow=annRow, annCol=annCol, main = main, info = info)	
	}
)
#' Plots a heatmap of the connectivity matrix of an KINOMO model.
setMethod('consensusmap', 'KINOMO', 
	function(object, ...){
		consensusmap(connectivity(object), ...)		
	}
)
#' Main method that redefines default values for arguments of \code{\link{aheatmap}}.
setMethod('consensusmap', 'matrix', 
	function(object, color='-RdYlBu'
			, distfun = function(x) as.dist(1-x), hclustfun = 'average'
			, Rowv = TRUE, Colv = "Rowv"
			, main = if( is.null(nr) || nr > 1 ) 'Consensus matrix' else 'Connectiviy matrix'
			, info = FALSE
			, ...){
				
		nr <- nrun(object)
		nb <- nbasis(object)
		info <- if( isTRUE(info) ){
					info <- NULL
					if( !is.null(nr) ) info <- c(info, paste("nrun:", nr))
					if( !is.null(nb) ) info <- c(info, paste("nbasis:", nb))
					info <- c(info, paste("cophcor:", round(cophcor(object), 3)))
				}else if( isFALSE(info) ) NULL
				else info
			
		aheatmap(object, color = color, ...
				, distfun = distfun, hclustfun = hclustfun
				, Rowv = Rowv, Colv = Colv
				, main = main
				, info = info)
	}
)

setOldClass('KINOMO.rank')
#' Draw a single plot with a heatmap of the consensus matrix obtained for each value of the rank, 
#' in the range tested with \code{\link{KINOMOEstimateRank}}.
#' 
#' @rdname KINOMO-compare
setMethod('consensusmap', 'KINOMO.rank', 
	function(object, ...){

		# plot the list of consensus matrix (set names to be used as default main titles)
		consensusmap(setNames(object$fit, paste("rank = ", lapply(object$fit, nbasis))), ...)
	}
)

#' 
#' @rdname KINOMO-compare
setMethod('consensusmap', 'list', 
	function(object, layout
			, Rowv = FALSE, main = names(object)
			, ...){
				
		opar <- par(no.readonly=TRUE)
		on.exit(par(opar))
		
		# define default layout
		if (missing(layout) ){
			n <- length(object)
			nr <- nc <- floor(sqrt(n))
			if( nr^2 != n ){
				nc <- nr + 1
				if( nr == 1 && nr*nc < n )
					nr <- nr + 1
			}
			
			layout <- c(nr, nc)		
		}
		if( !is.matrix(layout) ){
			if( !is.numeric(layout) )
				stop("invalid layout specification: must be a matrix or a numeric")
			if( length(layout) == 1 )
				layout <- c(layout, layout)
			layout <- matrix(1:(layout[1]*layout[2]), layout[1], byrow=TRUE)
		}
		
		graphics::layout(layout)
		res <- sapply(seq_along(object), function(i, ...){
			x <- object[[i]]
			
			# set main title
			main <- if( !is.null(main) && length(main) > 1 ){
				if( length(main) != length(object) )
					stop("consensusmap - Invalid length for argument `main`: should be either a single character string, or a list or vector of same length as ", deparse(substitute(object)))
				main[[i]]
			}			
			
			# call method for the fit
			consensusmap(x, ..., Rowv=Rowv, main=main)
		}, ...)
		invisible(res)
	}
)

#' Plots a heatmap of the basis matrix of the best fit in \code{object}.
setMethod('basismap', signature(object='KINOMOfitX'),
	function(object, ...){
		# call the method on the best fit
		basismap(minfit(object), ...)	
	}
)


setMethod('coefmap', signature(object='KINOMOfitX'),
	function(object
			, Colv=TRUE
			, annRow=NA, annCol=NA
			, tracks=c('basis', 'consensus:')
			, ...){
		
		x <- minfit(object)
		
		# process annotation tracks
		ptracks <- process_tracks(x, tracks, annRow, annCol)
		annRow <- ptracks$row 
		annCol <- ptracks$col
		# set special annotation handler
		specialAnnotation(2L, 'consensus', function() predict(object, what='consensus'))
		# row track handler is added in coefmap,KINOMO
		#
		
		## process ordering
		if( isString(Colv) ){
			if( Colv %in% c('consensus', 'cmap') )
				Colv <- consensushc(object, 'consensus')
		}
		##
		# call the method on the best fit
		coefmap(x, ..., Colv=Colv, annRow=annRow, annCol=annCol, tracks=NA)	
	}
)


setGeneric('cophcor', function(object, ...) standardGeneric('cophcor') )

#' 
setMethod('cophcor', signature(object='matrix'),
	function(object, linkage='average'){
		
		# check for empty matrix
		if( nrow(object)==0  || ncol(object)==0 )
		{
			warning("KINOMO::cophcor - NA produced [input matrix is of dimension ", nrow(object), "x", ncol(object), "]"
					, call.=FALSE)
			return(NA)
		}
		
		# safe-guard for diagonal matrix: to prevent error in 'cor'
		if( all(object[upper.tri(object)]==0) && all(diag(object)==object[1,1]) )
			return(1)
		# convert consensus matrix into dissimilarities
		d.consensus <- as.dist(1 - object)
		# compute cophenetic distance based on these dissimilarities
		hc <- hclust(d.consensus, method=linkage)
		d.coph <- cophenetic(hc)
		
		# return correlation between the two distances
		res <- cor(d.consensus, d.coph, method='pearson')
		return(res)
	}
)
#' Computes the cophenetic correlation coefficient on the consensus matrix
#' of \code{object}.
#' All arguments in \code{...} are passed to the method \code{cophcor,matrix}.
setMethod('cophcor', signature(object='KINOMOfitX'),
	function(object, ...){
		# compute the consensus matrix
		C <- consensus(object)
		
		return( cophcor(C, ...))
	}
)


