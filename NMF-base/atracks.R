setOldClass('annotationTrack')

#' Annotation Tracks
#' 
#' \code{.atrack} is an S4 generic method that converts an object into 
#' an annotation track object.
#' It provides a general and flexible annotation framework that is used 
#' by \code{\link{aheatmap}} to annotates heatmap rows and columns.
#' 
#' Methods for \code{.atrack} exist for common type of objects, which  
#' should provide enough options for new methods to define how annotation 
#' track are extracted from more complex objects, by coercing/filtering 
#' them into a supported type.
#' 
#' @param object an object from which is extracted annotation tracks
#' @param ... extra arguments to allow extensions and passed to the next method
#' call.
#' For \code{atrack}, arguments in \code{...} are concatenated into a single
#' \code{annotationTrack} object. 
#'
#' @rdname atrack 
#' @export 
#' @inline
#' @keywords internal
setGeneric('.atrack', function(object, ...) standardGeneric('.atrack'))

#' \code{is.atrack} tests if an object is an \code{annotationTrack} object.
#' 
#' @param x an R object
#' 
#' @rdname atrack
is.atrack <- function(x) is(x, 'annotationTrack')

aname <- function(x, name){
	return(x)
	if( missing(name) ){
		cn <- colnames(x)
		an <- attr(x, 'aname')
		name <- 
				if( !is.null(cn) ) cn		
				else if( !is.null(an) ) an
				else class(x)[1]
		attr(x, 'aname') <- name
	}else{
		attr(x, 'aname') <- name
		x
	}
}

#' \code{adata} get/sets the annotation parameters on an object
#' 
#' @param value replacement value for the complete annotation data list 
#' 
#' @rdname atrack
adata <- function(x, value, ...){
	if( missing(value) ){
		ad <- attr(x, 'annotationData')
		if( is.null(ad) ) ad <- list()
		
		# either return the annotationData itself or set values and return the object
		if( nargs() == 1L ) ad
		else{
			ad <- c(list(...), ad)
			ad <- ad[!duplicated(names(ad))]
			adata(x, ad)
		}
	}else{
		if( !is.list(value) )
			stop("Annotation data must be a list.")
		attr(x, 'annotationData') <- value
		x
	}
}

#' \code{amargin} get/sets the annotation margin, i.e. along which dimension of 
#' the data the annotations are to be considered.
#' 
#' @rdname atrack
amargin <- function(x, value){
	if( missing(value) ) adata(x)$margin
	else adata(x, margin=value)
}

#' \code{anames} returns the reference margin names for annotation tracks, 
#' from their embedded annotation data object.
#' 
#' @rdname atrack
anames <- function(x, default.margin){
	
	if( is.numeric(x) && length(x) == 1L ) NULL
	else if( is.vector(x) ) names(x)
	else{
		m <- amargin(x)
		if( is.null(m) && !missing(default.margin) ) m <- default.margin
		
		# special case for ExpressionSet  objects whose dimnames method returns NULL
		if( is(x, 'ExpressionSet') ) x <- Biobase::exprs(x)
		
		if( !is.null(m) ) dimnames(x)[[m]]
		else NULL
	}
}

#' \code{alength} returns the reference length for annotation tracks, 
#' from their embedded annotation data object
#' 
#' @param default.margin margin to use if no margin data is stored in the 
#' \code{x}.
#' 
#' @rdname atrack
alength <- function(x, default.margin){
	
	if( is.numeric(x) && length(x) == 1L ) as.integer(x)
	else if( is.vector(x) ) length(x)
	else{
		m <- amargin(x)
		if( is.null(m) && !missing(default.margin) ) m <- default.margin
		if( !is.null(m) ) dim(x)[m]
		else NULL
	}
	
}

test.match_atrack <- function(){
	
    checkEquals <- RUnit::checkEquals
	na <- paste("name_", 1:10, sep='') 
	mat <- as.matrix(setNames(1:10, na))
	
	.check <- function(x){		
		cat(class(x), " [", str_out(x, Inf, use.names=TRUE), "] :\n")
		y <- match_atrack(x, mat)
		print(y)
		checkEquals( class(y), class(x), "Same class as input")
		checkEquals( length(y), nrow(mat), "Correct length")
		checkEquals( names(y), rownames(mat), "Correct names")		
	}
	.test <- function(x){		
		
		.check(x)
		.check(sample(x))
		.check(x[1:5])
		.check(sample(x)[1:5])
		
		.check(setNames(x, na))
		.check(sample(setNames(x, na)))
		.check(setNames(x, rev(na)))		
		.check(setNames(x, na)[1:5])
		.check(setNames(x, na)[3:6])
		.check(setNames(x, na)[c(3,2,6)])
		
		x2 <- setNames(c(x[1:5], x[1:3]), c(na[1:5], paste("not_in_", 1:3, sep='')))
		.check(x2)
	}
	
	.test(letters[1:10])
	.test(1:10)
	.test(as.numeric(1:10) + 0.5)
	.test(c(rep(TRUE, 5), rep(FALSE, 5)))
	.test(factor(gl(2,5,labels=c("A", "B"))))
}

#' Extending Annotation Vectors
#' 
#' Extends a vector used as an annotation track to match the number of rows 
#' and the row names of a given data.
#' 
#' @param x annotation vector
#' @param data reference data
#' @return a vector of the same type as \code{x}
#' @export
#' 
match_atrack <- function(x, data=NULL){
	
	if( is.null(data) || length(x) == 0L ) return(x)
	
	# reorder and extend if a reference data matrix is provided
	refnames <- anames(data, default.margin=1L)
	reflength <- alength(data, default.margin=1L)
	
	# if no ref length (=> no refnames either): do nothing
	if( is.null(reflength) ) return(x)
	
	# special handling of character vectors
	if( is.character(x) && is.null(names(x)) && !is.null(refnames) ){
#		if( !any(names(x) %in% refnames) && any(x %in% refnames) ){
		if( any(x %in% refnames) ){
			vmessage("match_atrack - Annotation track [", str_out(x, 3, use.names=TRUE), "] has some values matching data names: converting into a logical using values as names.")
			x <- setNames(rep(TRUE, length(x)), x)
		}
	}
	
	# reorder based on names
	.hasNames <- FALSE
	if( !is.null(names(x)) && !is.null(refnames) ){
		inref <- names(x) %in% refnames
		if( !all(inref) ){
			vmessage("match_atrack - Annotation track [", str_out(x, 3, use.names=TRUE), "] has partially matching names: subsetting track to match data")
			x <- x[inref]
			if( length(x) == 0L )
				vmessage("match_atrack - Subset annotation track is empty")
		}else
			vmessage("match_atrack - Annotation track [", str_out(x, 3, use.names=TRUE), "] using names as identifiers")
		.hasNames <- TRUE
					
		if( anyDuplicated(names(x)) ){
			dups <- duplicated(names(x))
			vmessage("match_atrack - Annotation track [", str_out(x, 3, use.names=TRUE), "]: removing duplicated names [", str_out(x[dups], 3, use.names=TRUE),"]")
			x <- x[!dups]
		}
	}
	
	lx <- length(x)
	if( lx > reflength ){
		stop("match_atrack - Invalid annotation track [", str_out(x, 3, use.names=TRUE), "]: more elements [", lx, "] than rows in data [", reflength, "].")
	}
	if( lx == reflength ){
		# reorder if necessary
		res <- 
			if( !.hasNames ) x
			else x[match(refnames, names(x))]
		return(res)
	}
	
	# build similar vector of correct size
	res <- 
		if( is.factor(x) ) setNames(factor(c(x, rep(NA, reflength-lx)), levels=c(levels(x), NA)), refnames)
		else setNames(c(x, rep(NA, reflength-lx)), refnames)
	res[1:lx] <- NA
	
	# if not using names
	if( !.hasNames ){						
		if( is.integer(x) ) res[x] <- x		
		else res[1:lx] <- x
	}else{
		# put the values of x at the write place
		res[match(names(x), refnames)] <- x
	}
	res
}

#' The default method converts character or integer vectors into factors. 
#' Numeric vectors, factors, a single NA or \code{annotationTrack} objects are 
#' returned unchanged (except from reordering by argument \code{order}).
#' Data frames are not changed either, but class 'annotationTrack' is appended 
#' to their original class set.
#' 
#' @param data object used to extend the annotation track within a given data 
#' context.
#' It is typically a matrix-like object, against which annotation specifications 
#' are matched using \code{\link{match_atrack}}.
#'  
#' 
setMethod('.atrack', signature(object='ANY'),
	function(object, data=NULL, ...){
		
		# recursive on list
		if( is.list(object) ){
			object <- object[!sapply(object, function(x) length(x) == 0 || is_NA(x) )]
			res <- 
					if( length(object) == 0 ) NULL
					else{
						# convert into a list of tracks
						sapply(object, .atrack, data=data, ..., simplify=FALSE) 
					}
			return(res)
			
		}else if( is.null(object) || is_NA(object) || is.atrack(object) ) object
		else{
			# extend to match the data
			object <- match_atrack(object, data)
			
			# apply convertion rules for standard classes
			if( is.logical(object) ) aname(as.factor(ifelse(object, 1, NA)), "Flag")
			else if( is.integer(object) ){
				if( any(wna <- is.na(object)) )
					aname(as.factor(ifelse(!wna, 1,NA)), "Flag")
				else 
					aname(as.numeric(object), "Level")
			} 
			else if( is.character(object) ) aname(as.factor(object), "Group")
			else if( is.factor(object) ) aname(object, "Factor")
			else if( is.numeric(object) ) aname(object, "Variable")				
			else stop("atrack - Invalid annotation item `"
						, substitute(object)
						, "`: must be a factor, or a logical, character, numeric or integer vector")
		}
		
	}
)

#' Adds annotation tracks resolving special track codes that are those elements that start with ":"
setMethod('.atrack', 'character', 
	function(object, ...){
		
		# check for special escaped track code
    	if( length(i <- atrack_code(object)) ){
			if( length(object) == 1L ) object
			else if( length(i) == length(object) ) as.list(object)
			else{
#				spe <- object[i]
#				object <- sub("^\\\\:", ":", object[-i])
#				t <- callNextMethod() 
#				c(list(t), spe)
				callNextMethod()
			}
		}else{
#			object <- sub("^\\\\:", ":", object)
			callNextMethod()
		}
	}
)
#' Shorcut for `.atrack(as.data.frame(object), ...)`
setMethod('.atrack', 'matrix', function(object, ...) .atrack(as.data.frame(object), ...) )
#' Considers each column as an annotation variable.
#' Equivalent to `.atrack(as.list(object), ...)`
setMethod('.atrack', 'data.frame', function(object, ...) .atrack(as.list(object), ...) )

# tells if an object is a special annotation track code
is_track_code <- function(x) isString(x) && grepl("^[:$]", x)

atrack_code <- function(x, value=FALSE){
	
	# check each track item
	ac <- sapply(x, is_track_code)
	i <- which(ac)
	
	if( !value ) i # return indexes
	else if( length(i) ) unlist(x[i]) # return values
	
}


match_atrack_code <- function(x, table, ...){
	# pre-pend ':'
	table.plain <- sub("^:", '', table)	
	table <- str_c(':', table.plain)
	
	# convert into an annotation track
	if( !is.atrack(x) ) x <- atrack(x, ...)
	
	m <- sapply(x, function(x){
		if( isString(x) ) charmatch(x, table, nomatch=0L)
		else 0L
	})
	
	if( length(i <- which(m!=0L)) ){
		if( is.null(names(m)) ) names(m) <- rep('', length(m))
		names(m)[i] <- table.plain[m[i]]
	}
	m
}

#' \code{atrack} creates/concatenates \code{annotationTrack} objects
#' 
#' @param order an integer vector that indicates the order of the annotation
#' tracks in the result list
#' @param enforceNames logical that indicates if missing track names should 
#' be generated as \code{X<i>}
#' @param .SPECIAL an optional list of functions (with no arguments) that are 
#' called to generate special annotation tracks defined by codes of the form 
#' \code{':NAME'}.
#' e.g., the function \code{link{consensusmap}} defines special tracks 
#' \code{':basis'} and \code{':consensus'}.
#' 
#' If \code{.SPECIAL=FALSE}, then any special tracks is discarded and a warning
#' is thrown. 
#'    
#' @param .DATA data used to match and extend annotation specifications.
#' It is passed to argument \code{data} of the \code{.atrack} methods, which 
#' in turn use pass it to \code{\link{match_atrack}}.
#' 
#' @param .CACHE an \code{annotationTrack} object with which the generated
#' annotation track should be consistent.
#' This argument is more for internal/advanced usage and should not be used 
#' by the end-user.  
#' 
#' @return \code{atrack} returns a list, decorated with class 
#' \code{'annotationTrack'}, where each element contains the description 
#' of an annotation track.
#'  
#' @rdname atrack
#' @export 
atrack <- function(..., order = NULL, enforceNames=FALSE, .SPECIAL=NA, .DATA = NULL, .CACHE = NULL){
	
	# cbind object with the other arguments
	l <- list(...)
	if( length(l) == 1L && is.atrack(l[[1]]) )
		object <- l[[1L]]
	else if( length(l) > 0 ){
		
		object <- list()
		#print(l)
		lapply(seq_along(l), function(i){
					x <- l[[i]]
					if( is_NA(x) || is.null(x) )
						return()
					
					xa <- .atrack(x, data=.DATA)
					
					if( is_NA(xa) || is.null(xa) )
						return()
					
					n <- names(object)
					# convert into a list
					if( !is.list(xa) )
						xa <- setNames(list(xa), names(l)[i])
							
					# remove NA and NULL elements
					if( is.null(xa) || is_NA(xa) ) return()
					# cbind with previous tracks
					if( is.null(object) ) object <<- xa
					else object <<- c(object, xa)
					
				})
	}
	
	# exit now if object is NULL
	if( is.null(object) ) return()
	if( !length(object) ) return( annotationTrack() )
	
	# add class 'annotationTrack' if not already there 
	# (needed before calling match_atrack_code)
	object <- annotationTrack(object)
	
	# substitute special tracks
	if( is.list(.SPECIAL) ){
	
#		str(object)
		m <- match_atrack_code(object, names(.SPECIAL))
		i_spe <- which(m!=0L)
		if( length(i_spe) ){
			# add names where needed
			if( is.null(names(object)) ) names(object) <- rep('', length(object))
			
			# remove duplicated special tracks
			if( anyDuplicated(m[i_spe]) ){
				# enforce name consistency if necessary
				g <- split(i_spe, m[i_spe])
				sapply(g, function(i){
					n <- names(object)[i]
					if( length(n <- n[n!='']) )					
						names(object)[i] <<- n[1L] 
				})
				#
				idup <- which(duplicated(m) & m!=0L)
				object <- object[-idup]
				m <- m[-idup]
				i_spe <- which(m!=0L)
			}
			#
		
			# enforce names consistent with the CACHE
			if( anyValue(.CACHE) ){
				if( !is.atrack(.CACHE) )
					stop("Argument .CACHE should be an annotation track object. [", class(.CACHE), ']')
				i_spe_cache <- atrack_code(.CACHE)
				if( length(i_spe_cache) ){
					.CACHE_SPE <- unlist(.CACHE[i_spe_cache])
					if( !is.null(names(.CACHE_SPE)) ){
						sapply(i_spe, function(i){
							x <- object[[i]]
							if( names(object)[i] == '' 
								&& !is_NA(j <- match(x, .CACHE_SPE)) 
								&& names(.CACHE_SPE)[j] != ''){
								names(object)[i] <<- names(.CACHE_SPE)[j]
							}
						})
					}
				}
			}
			# compute value
			a <- sapply(m[i_spe], function(i) .SPECIAL[[i]](), simplify=FALSE)
			object[i_spe] <- a # NB: this does not change the names
			# reset names
			nm <- names(object)[i_spe]
			names(object)[i_spe] <- ifelse(nm!='', nm, names(a))
		}
		
		# remove special tracks if necessary
		if( length(i <- atrack_code(object)) ){
			warning("Discarding unresolved special annotation tracks: "
					, str_out(unlist(object[i]), use.names=TRUE))
			object <- object[-i] 
		}
	}
	
	# generate names
	if( enforceNames ){
		n <- names(object)
		xnames <- paste('X', 1:length(object), sep='')
		if( is.null(n) ) names(object) <- xnames
		else names(object)[n==''] <- xnames[n==''] 
	}
	
	# reorder if necessary
	if( !is.null(order) ){
		object <- sapply(object, function(x) x[order], simplify=FALSE)
		#lapply(seq_along(object), function(i) object[[i]] <<- object[[i]][order])
	}
	
	#print(object)
	# return object
	annotationTrack(object)
}

#' \code{annotationTrack} is constructor function for \code{annotationTrack} object
#' 
#' @rdname atrack
annotationTrack <- function(x = list()){
	if( !is.atrack(x) )
		class(x) <- c('annotationTrack', if( nargs() ) class(x))
	x
} 

#setGeneric('atrack<-', function(object, value) standardGeneric('atrack<-'))
#setReplaceMethod('atrack', signature(object='ANY', value='ANY'),
#	function(object, value){
#		# if the annotation track is not NA: convert it into a atrack 
#		# and set the value
#		if( !is_NA(object) && length(value) > 0 ){
#			object <- atrack(object, value)
#		}
#		object
#	}
#)
#setReplaceMethod('atrack', signature(object='annotationTrack'),
#	function(object, value, replace = FALSE){
#		if( !replace && length(value) > 0 )	atrack(object, value)			
#		else if( replace ) atrack(value)
#		else object
#	}
#)
