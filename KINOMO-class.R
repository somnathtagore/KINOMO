#library(R.utils)

#' @include utils.R
#' @include versions.R
#' @include algorithmic.R
#' @include aheatmap.R
NULL


#' @rdname advanced
#' @name advanced-KINOMO
NULL 

# declare old S3 class 'proc_time' to use it as a slot for class KINOMO 
setOldClass('proc_time', prototype=numeric())

################################
# Class: KINOMO
################################






#' 
setClass('KINOMO'
		, representation(
			misc = 'list' # misceleneaous data used during fitting
		)
		, contains = 'VIRTUAL')



#' @export
setGeneric('fitted', package='stats')
 
setMethod('fitted', signature(object='KINOMO'),
		function(object, ...){
			stop("KINOMO::fitted is a pure virtual method of interface 'KINOMO'. It should be overloaded in class '", class(object),"'.")
		}
)


#'  
setGeneric('basis', function(object, ...) standardGeneric('basis') )
#' Default method returns the value of S3 slot or attribute \code{'basis'}.
#' It returns \code{NULL} if none of these are set. 
#' 
#' Arguments \code{...} are not used by this method.
setMethod('basis', signature(object='ANY'),
	function(object, ...){
		if( is.list(object) && 'basis' %in% names(object) ) object[['basis']]
		else attr(object, 'basis')
	}
)
#' @param all a logical that indicates whether the complete matrix factor 
#' should be returned (\code{TRUE}) or only the non-fixed part.
#' This is relevant only for formula-based KINOMO models that include fixed basis or 
#' coefficient terms.
#' 
setMethod('basis', signature(object='KINOMO'),
	function(object, all=TRUE, ...){
		if( all || !length(i <- ibterms(object)) ){
			# return all coefficients
			.basis(object, ...)
		} else {
			# remove fixed basis
			.basis(object, ...)[, -i]
		}
	}
)



#' @export
setGeneric('.basis', function(object, ...) standardGeneric('.basis') )

setMethod('.basis', signature(object='KINOMO'),
	function(object, ...){
		stop("KINOMO::.basis is a pure virtual method of interface 'KINOMO'. It should be overloaded in class '", class(object),"'.")
	}
)

#' @export
#' @rdname basis-coef-methods

setGeneric('basis<-', function(object, ..., value) standardGeneric('basis<-') )

#' 
setReplaceMethod('basis', signature(object='KINOMO', value='ANY'), 
	function(object, use.dimnames = TRUE, ..., value){
		
        # error if passed extra arguments
        if( length(xargs<- list(...)) ){
            stop("basis<-,KINOMO - Unused arguments: ", str_out(xargs, Inf, use.names = TRUE))
        }
        
        # backup old dimnames to reapply them on exit
        if( !use.dimnames ) odn <- dimnames(object)
        nb_old <- nbasis(object)
        
		# only set non-fixed terms
		if( !nbterms(object) ) .basis(object) <- value
		else{
			i <- ibasis(object)
			.basis(object)[,i] <- value[, i]
		}
        # adapt coef if empty
        if( !hasCoef(object) ){
            x <- basis(object)
            .coef(object) <- rbind(coef(object)[1:min(nb_old, ncol(x)), , drop = FALSE], matrix(NA, max(ncol(x)-nb_old, 0), 0))
#            .coef(object) <- coef(object)[1:ncol(x), , drop = FALSE] 
        }
		# check object validity
		validObject(object)
        
        # update other factor if necessary
        if( use.dimnames ) basisnames(object) <- colnames(basis(object))
        else if( !length(odn) ) dimnames(object) <- NULL
        else dimnames(object) <- mapply(head, odn, dim(object), SIMPLIFY = FALSE)
        
		object
	}	
)
#' @param value replacement value 
#' @rdname basis-coef-methods
#' @export
setGeneric('.basis<-', function(object, value) standardGeneric('.basis<-') )

setReplaceMethod('.basis', signature(object='KINOMO', value='matrix'), 
	function(object, value){ 
		stop("KINOMO::.basis<- is a pure virtual method of interface 'KINOMO'. It should be overloaded in class '", class(object),"'.")
	} 
)

#' @export
setGeneric('loadings', package='stats')


#' 
#' @rdname basis-coef-methods
setMethod('loadings', 'KINOMO', function(x) basis(x) )


#' @export
setGeneric('coef', package='stats')

setMethod('coef', 'KINOMO',
	function(object, all=TRUE, ...){
		
		if( all || !length(i <- icterms(object)) ){
			# return all coefficients
			.coef(object, ...)
		} else {
			# remove fixed coefficients
			.coef(object, ...)[-i, ]
		}
		
	}
)


#' @export
setGeneric('.coef', function(object, ...) standardGeneric('.coef'))

setMethod('.coef', signature(object='KINOMO'),
	function(object, ...){
		stop("KINOMO::.coef is a pure virtual method of interface 'KINOMO'. It should be overloaded in class '", class(object),"'.")
	}
)

#' @export
#' @rdname basis-coef-methods

setGeneric('coef<-', function(object, ..., value) standardGeneric('coef<-') )
#' Default methods that calls \code{.coef<-} and check the validity of the 
#' updated object. 
setReplaceMethod('coef', signature(object='KINOMO', value='ANY'), 
	function(object, use.dimnames = TRUE, ..., value){
		
        # error if passed extra arguments
        if( length(xargs<- list(...)) ){
            stop("coef<-,KINOMO - Unused arguments: ", str_out(xargs, Inf, use.names = TRUE))
        }
        # backup old dimnames to reapply them on exit
        if( !use.dimnames ) odn <- dimnames(object)
        nb_old <- nbasis(object)
        
		# only set non-fixed terms
		if( !ncterms(object) ) .coef(object) <- value
		else{
			i <- icoef(object)
			.coef(object)[i, ] <- value[i, ]
		}
        # adapt basis if empty before validation
        if( !hasBasis(object) ){
            x <- coef(object)
            .basis(object) <- cbind(basis(object)[, 1:min(nb_old, nrow(x)), drop = FALSE], matrix(NA, 0, max(nrow(x)-nb_old, 0))) 
        }
		# check object validity
		validObject(object)
        
        # update other factor if necessary
        if( use.dimnames ) basisnames(object) <- rownames(coef(object))
        else if( !length(odn) ) dimnames(object) <- NULL
        else dimnames(object) <- mapply(head, odn, dim(object), SIMPLIFY = FALSE)
        
            
		object
	}	
)

#' @export
#' @rdname basis-coef-methods
setGeneric('.coef<-', function(object, value) standardGeneric('.coef<-') )

setReplaceMethod('.coef', signature(object='KINOMO', value='matrix'), 
	function(object, value){ 
		stop("KINOMO::.coef<- is a pure virtual method of interface 'KINOMO'. It should be overloaded in class '", class(object),"'.")
	} 
)


#' @rdname basis-coef-methods
setGeneric('coefficients', package='stats')
#' Alias to \code{coef,KINOMO}, therefore also pure virtual.

setMethod('coefficients', signature(object='KINOMO'), selectMethod('coef', 'KINOMO'))


#' 
setGeneric('scoef', function(object, ...) standardGeneric('scoef') )

setMethod('scoef', 'KINOMO',
	function(object, scale=1){
		sweep(coef(object), 2L, colSums(coef(object)) / scale, '/')
	}
)

setMethod('scoef', 'matrix',
	function(object, scale=1){
		sweep(object, 2L, colSums(object) / scale, '/')
	}
)

unit.test(scoef, {
	x <- rKINOMO(3, 10, 5)
	checkIdentical(colSums(scoef(x)), rep(1, nbasis(x))
		, "Default call: columns are scaled to sum-up to one")
	checkIdentical(colSums(scoef(x, 100)), rep(1, nbasis(x))
		, "Scale=10: columns are scaled to sum-up to 10")
})


#' 
scale.KINOMO <- function(x, center=c('basis', 'coef'), scale=1){
	
	# determine base value
	if( missing(center) ) center <- match.arg(center)
	base <- center
	delta <-
	if( is.character(base) ){
		base <- match.arg(center)
		if( base == 'basis' ) colSums(basis(x))
		else{
			scale <- 1/scale 
			1 / rowSums(coef(x))
		}	
	}else if( is.numeric(base) ) base
	else stop("Invalid base value: should be a numeric or one of "
			, str_out(c('none', 'basis', 'coef')))

	# scale
	D <- scale/delta 
	# W <- W * D
	basis(x) <- sweep(basis(x), 2L, D, '*')
	# H <- D^-1 * H
	coef(x) <- sweep(coef(x), 1L, D, '/')
	x
	
}

unit.test("scale", {
			
	r <- 3
	x <- rKINOMO(r, 10, 5)
	
	.lcheck <- function(msg, rx, ref, target){
		.msg <- function(...) paste(msg, ':', ...)
		checkTrue(!identical(basis(x), basis(rx)), .msg("changes basis matrix"))
		checkTrue(!identical(coef(x), coef(rx)), .msg("changes coef matrix"))
		checkEqualsNumeric(fitted(x), fitted(rx), .msg("fitted target is identical"))
		
		brx <- colSums(basis(rx))
		crx <- rowSums(coef(rx)) 
		if( target == 1 ){
			checkEquals(brx, ref, .msg("correctly scales basis components"))
			checkTrue(!all(crx==ref), .msg("does not force scale on coefficient matrix"))
		}else{
			checkTrue(!all(brx==ref), .msg("does not force scale on basis matrix"))
			checkEquals(crx, ref
					, .msg("correctly scales rows of coef matrix"))
		}
	}
	
	.check <- function(msg, ref, ...){
		.lcheck(str_c(msg, " + argument center='basis'")
				, scale(x, center='basis', ...), ref, 1)
		.lcheck(str_c(msg, " + argument center='coef'")
				, scale(x, center='coef', ...), ref, 2)
	}
	
	.lcheck("Default call", scale(x), rep(1, r), 1)
	.check("Missing argument scale", rep(1, r))
	.check("Argument scale=10", rep(10, r), scale=10)
	s <- runif(r)
	.check("Argument scale=numeric", s, scale=s)
 
})



#' @family KINOMO-interface
setGeneric('rKINOMO', function(x, target, ...) standardGeneric('rKINOMO') )

# Define the loading namespace
.PKG.NAMESPACE <- packageEnv()


#' 
is.KINOMO <- function(x){
	
	# load definition for base class KINOMO
	clref <- getClass('KINOMO', .Force=TRUE, where=.PKG.NAMESPACE)
	is(x, clref)
	
}

unit.test(is.KINOMO,{
	checkTrue(!is.KINOMO(1:4), "on vector: FALSE")
	checkTrue(!is.KINOMO(list(1:4)), "on list: FALSE")
	checkTrue(is.KINOMO('KINOMO'), "on 'KINOMO': TRUE")
	checkTrue(is.KINOMO('KINOMOstd'), "on 'KINOMOstd': TRUE")
	checkTrue( is.KINOMO( KINOMOModel(3) ), "on empty model: TRUE")
	checkTrue( is.KINOMO( rKINOMO(3, 20, 10) ), "on random model: TRUE")
	checkTrue( is.KINOMO( KINOMO(rmatrix(20,10), 3) ), "on KINOMOfit object: TRUE") 
})

isKINOMOclass <- function(x){
	
	if( is.character(x) ){ # test object is a class that extends KINOMO
		# load definition for base class KINOMO
		clref <- getClass('KINOMO', .Force=TRUE, where=.PKG.NAMESPACE)
		cl <- getClass(x, .Force=TRUE, where=.PKG.NAMESPACE)
		if( is.null(cl) )
			cl <- getClass(x, .Force=TRUE)
		extends(cl, clref)
	}else 
		FALSE
	
}

################################

# Taken from Biobase
selectSome <- function (obj, maxToShow = 5) 
{
    len <- length(obj)
    if (maxToShow < 3) 
        maxToShow <- 3
    if (len > maxToShow) {
        maxToShow <- maxToShow - 1
        bot <- ceiling(maxToShow/2)
        top <- len - (maxToShow - bot - 1)
        nms <- obj[c(1:bot, top:len)]
        c(as.character(nms[1:bot]), "...", as.character(nms[-c(1:bot)]))
    }
    else if (is.factor(obj)) 
        as.character(obj)
    else obj
}


.showFixedTerms <- function(x, ...){
	
	s <- 
	sapply(x, function(t){
		s <- 
		if( is.factor(t) ) selectSome(levels(t), ...)
		else selectSome(t, ...)
		s <- str_out(s, Inf, quote=FALSE)
		if( is.factor(t) )
    		s <- str_c('<', s, ">")
		s
	})
	paste(names(s), '=', s)
}

#' Show method for objects of class \code{KINOMO}
#' @export
setMethod('show', 'KINOMO', 
		function(object)
		{
			cat("<Object of class:", class(object), ">\n", sep='')
			cat("features:", nrow(object), "\n")
			cat("basis/rank:", nbasis(object), "\n")
			cat("samples:", ncol(object), "\n")
			# show fixed terms
			if( (n <- ncterms(object)) ){
				cat("fixed coef [", n, "]:\n" 
					, str_c('  ', .showFixedTerms(cterms(object), 4), collapse="\n")
					, "\n", sep='')
			}
			if( (n <- nbterms(object)) ){
				cat("fixed basis [", n, "]:\n" 
						, str_c('  ', .showFixedTerms(bterms(object), 4), collapse="\n")
						, "\n", sep='')
			}
			# show the miscellaneous model parameters
			if( length(object@misc) > 0L ){
				cat("miscellaneous:", str_desc(object@misc, exdent=12L), ". (use 'misc(object)')\n")
			}
		}
)



#'    
setMethod('nbasis', signature(x='ANY'), 
	function(x, ...)
	{
		if( !is.null(n <- ncol(basis(x, ...))) ) n
		else if( is.list(x) && 'nbasis' %in% names(x) ) x[['nbasis']]
		else attr(x, 'nbasis')
	}
)


#' @export
setMethod('dim', signature(x='KINOMO'), 
	function(x){
		c(nrow(basis(x)), ncol(coef(x)), nbasis(x))	
	}
)




#' 
setGeneric('basisnames', function(x, ...) standardGeneric('basisnames') )

#' @rdname dimnames
setMethod('basisnames', signature(x='ANY'), 
	function(x)
	{
		colnames(basis(x))
	}
)

#' @export
#' @rdname dimnames 
setGeneric('basisnames<-', function(x, ..., value) standardGeneric('basisnames<-') )

#' @rdname dimnames 
setReplaceMethod('basisnames', 'ANY', 
	function(x, ..., value)
	{
		rownames(.coef(x)) <- value
		colnames(.basis(x)) <- value
		x
	}
)


#' @export
setMethod('dimnames', 'KINOMO', 
	function(x){
		b <- dimnames(basis(x))
		if( is.null(b) )
			b <- list(NULL, NULL)
		c <- dimnames(coef(x))
		if( is.null(c) )
			c <- list(NULL, NULL)
		l <- c(b[1],c[2],b[2])
		if( all(sapply(l, is.null)) ) NULL else l
	}
)

#' @rdname dimnames
#' @export
setReplaceMethod('dimnames', 'KINOMO', 
	function(x, value){
		if( !is.list(value) && !is.null(value) )
			stop("KINOMO::dimnames - Invalid value: must be a list or NULL.")
		
		if( length(value) == 0 )
			value <- NULL
		else if( length(value) == 1 )
			value <- c(value, list(NULL, NULL))			
		else if( length(value) == 2 ) # if only the two first dimensions reset the third one
			value <- c(value, list(NULL))
		else if( length(value)!=3 ) # check length of value
			stop("KINOMO::dimnames - invalid argument 'value' [a 2 or 3-length list is expected]")
		
		# only set relevant dimensions
		if( length(w <- which(dim(x) == 0)) ){
			value[w] <- sapply(value[w], function(x) NULL, simplify=FALSE)
		}
		# set dimnames 
		dimnames(.basis(x)) <- value[c(1,3)]
		dimnames(.coef(x)) <- value[c(3,2)]
		# return updated model
		x	
	}
)




#' 
setMethod('[', 'KINOMO', 
	function (x, i, j, ..., drop = FALSE)
	{
		k <- NULL
		mdrop <- missing(drop)
		
		# compute number of arguments: x and drop are always passed
		Nargs <- nargs() - !mdrop
		single.arg <- FALSE
		k.notmissing <- FALSE
		if( !missing(i) && Nargs < 3L ){
			k <- i
			single.arg <- TRUE
		}
		else if( Nargs > 3L ){
			dots <- list(...)
			if( length(dots) != 1 )
				warning("KINOMO::[ - using only the first extra subset index, the remaining ", length(dots)-1," are discarded.")
			k <- dots[[1]]
			k.notmissing <- TRUE
		}
		
		# no indice was provided => return the object unchanged
		if ( missing(i) && missing(j) && !k.notmissing ) {
			# check if there is other arguments
			if (length(list(...)) != 0)
				stop("KINOMO::[] method - please specify which features, samples or basis to subset. See class?KINOMO.")
			# otherwise return the untouched object
			return(x)
		}
		
		# subset the rows of the basis matrix		
		if ( !missing(i) && !single.arg )			
			.basis(x) <- basis(x)[i, , drop = FALSE]
		
		# subset the columns of mixture coefficient matrix		
		if (!missing(j)) 
			.coef(x) <- coef(x)[, j, drop = FALSE]
		
		# subset the basis: columns of basis matrix and row of mixture coefficient matrix		
		if( single.arg || k.notmissing ){
			.basis(x) <- basis(x)[, k, drop = FALSE]
			# return basis only single arg and drop=TRUE 
			if( single.arg && ((mdrop && length(k) == 1L) || drop) ) return( drop(basis(x)) )
			.coef(x) <- coef(x)[k, , drop = FALSE]
		}
		
		# if drop is TRUE and only one dimension is missing then return affected matrix
		if( !single.arg && drop ){
			if( missing(i) && !missing(j) )
				return( drop(coef(x)) )
			else if( missing(j) && !missing(i) )
				return( drop(basis(x)) )
		}
		
		# return subset object
		return(x)
	}
)


#' @export
misc <- function(object, ...){
	if( !isS4(object) && is.list(object) ) object[['misc']]
	else attr(object, 'misc')
}

#' @export 
setMethod('$', 'KINOMO', 
		function(x, name){ 
			x@misc[[name, exact=TRUE]]; 
		} 
)

#' @export
setReplaceMethod('$', 'KINOMO',
	function(x, name, value) {
		x@misc[[name]] <- value
		x
	}
)

#' @importFrom utils .DollarNames
setGeneric('.DollarNames', package='utils')

#' @method .DollarNames KINOMO
#' @export
.DollarNames.KINOMO <- function(x, pattern = "") grep(pattern, names(misc(x)), value=TRUE)

#' Auto-completion for \code{\linkS4class{KINOMO}} objects
#' @rdname KINOMO-class
#' @export
setMethod('.DollarNames', 'KINOMO', .DollarNames.KINOMO)


is.empty.KINOMO <- function(x, ...){
	nrow(x) == 0 && ncol(x) == 0
}


#' \code{hasBasis} tests whether an objects contains a basis matrix -- returned by 
#' a suitable method \code{basis} -- with at least one row.
#' 
#' @rdname types
#' @export
hasBasis <- function(x) nbasis(x) && nrow(basis(x)) != 0L

#' \code{hasBasis} tests whether an objects contains a coefficient matrix 
#' -- returned by a suitable method \code{coef} -- with at least one column.
#' 
#' @rdname types
#' @export
hasCoef <- function(x) nbasis(x) && ncol(coef(x)) != 0L

#' \code{is.partial.KINOMO} tests whether an KINOMO model object contains either an empty 
#' basis or coefficient matrix.
#' It is a shorcut for \code{!hasCoef(x) || !hasBasis(x)}.
#' 
#' @rdname types
#' @export
is.partial.KINOMO <- function(x) !hasCoef(x) || !hasBasis(x)


#' 
setMethod('rmatrix', 'KINOMO', 
	function(x, ...){
		a <- fitted(x)
		a + rmatrix(a, ...)
	}
)

unit.test('rmatrix,KINOMO',{
	
	x <- KINOMOModel(3, 20, 5)
	checTrue(is.matrix(y <- rmatrix(x)), "default call: no error")
	checkIdentical(dim(y), dim(x)[1:2], "default call: correct dimension")
	checkTrue( !any(is.na(basis(y))), 'default call: no NAs in basis anymore')
	checkTrue( !any(is.na(coef(y))), 'default call: no NAs in coef anymore')
	checkTrue( max( max(abs(basis(y)-basis(x))), max(abs(coef(y)-coef(x))) ) <= 1
			, "default call: max difference is <= 1")
	
	set.seed(123)
	y <- rmatrix(x)
	set.seed(123)
	ref <- matrix(runif(nrow(x)*ncol(x)), nrow(x))
	checkIdentical(ref, y - fitted(x), "default call: add uniform random noise to fitted matrix")
	
	set.seed(123)
	ref <- matrix(rnorm(nrow(x)*ncol(x)), nrow(x))
	set.seed(123)
	y <- rmatrix(x, rnorm)	
	checkIdentical(ref, y - fitted(x), "dist is taken into account: add normal random noise to fitted matrix")
	set.seed(123)
	y <- rmatrix(x, dist=rnorm)
	checkIdentical(ref, y - fitted(x), "dist is taken into account: add normal random noise to fitted matrix")
	
	set.seed(123)
	checTrue(is.matrix(y <- rmatrix(x, max=10)), "call with arg max=10: no error")
	checkTrue( max( max(abs(basis(y)-basis(x))), max(abs(coef(y)-coef(x))) ) <= 10
			, "call with arg max=10: max difference is 10")
	checkTrue( max( max(abs(basis(y)-basis(x))), max(abs(coef(y)-coef(x))) ) >= 5
			, "call with arg max=10: max difference is >= 5")
	
})



###% Utility function used to sets default elements in a list if they are
###% not already set
###% The default values are given in argument ...
.set.list.defaults <- function(input.list, ...){
	expand_list(input.list, ..., .exact=FALSE)
}

###% Partially match arguments for a given function 
.match.call.args <- function(x, fun, in.fun=NULL, call=NULL){
	stopifnot( is.character(fun) && length(fun) == 1 )
	if( length(x) == 0 ) return(x)
	x.ind <- charmatch(x, args <- formalArgs(getFunction(fun)))
	
	sapply(seq(length(x)), function(i){
			ind <- x.ind[i]
			# the argument is not part of the call: keep it unchanged
			if( is.na(ind) ) return(x[i])
			# multiple matches: error
			if( ind == 0 ){
				alt <- paste(grep(paste('^', x[i], sep=''), args, value=TRUE), collapse=', ')
				stop(if( !is.null(call) ) c(call, ' - '), "Multiple match for argument '", x[i], "' of function '"
					, if( is.null(in.fun) ) fun else in.fun, "' [use one of: ", alt, "]"
					, call.=FALSE)
			}
			# return the matched full names
			args[ind]
	})
}




#' 
setGeneric('summary', package='base')


setMethod('summary', signature(object='KINOMO'), 
		function(object, class, target){
			
			res <- numeric()
			
			## IMPORTANT: if adding a summary measure also add it in the sorting 
			## schema of method KINOMOList::summary to allow ordering on it
			
			# rank
			res <- c(res, rank=nbasis(object))
			# compute sparseness
			res <- c(res, sparseness=sparseness(object))
			
			# if class is provided: also computes entropy and purity
			if( !missing(class) ){
				# compute purity
				res <- c(res, purity=purity(object, class))
				# compute entropy
				res <- c(res, entropy=entropy(object, class))
			}
			
			# if the target is provided compute the RSS
			if( !missing(target) ){
				RSS <- rss(object, target)
				res <- c(res, rss=RSS)
				# explained variance
				res <- c(res, evar=evar(object, target))
			}
            
            # compute mean silhouette width
            siS <- silhouette(object, what = 'samples')
            siF <- silhouette(object, what = 'features')
            res <- c(res, silhouette.coef = if( !is_NA(siS) ) summary(siS)$avg.width else NA
                    , silhouette.basis = if( !is_NA(siF) ) summary(siF)$avg.width else NA)
			
			# return result
			return(res)
		}
)



setGeneric('sparseness', function(x, ...) standardGeneric('sparseness') )
#' Base method that computes the sparseness of a numeric vector.
#' 
#' It returns a single numeric value, computed following the definition 
#' given in section \emph{Description}. 
setMethod('sparseness', signature(x='numeric'), 
	function(x){
		# get length of x
		n <- length(x)
		# compute and return the sparseness
		( sqrt(n) - sum(abs(x)) / sqrt(sum(x^2)) ) / (sqrt(n)-1)
	}
)
#' Computes the sparseness of a matrix as the mean sparseness of its column vectors.
#' It returns a single numeric value.
setMethod('sparseness', signature(x='matrix'), 
	function(x){
		# compute the sparseness of each column
		s <- apply(x, 2, sparseness)
		
		# return the mean sparseness
		mean(s)
	}
)
#' Compute the sparseness of an object of class \code{KINOMO}, as the sparseness of 
#' the basis and coefficient matrices computed separately.
#' 
#' It returns the two values in a numeric vector with names \sQuote{basis} and \sQuote{coef}. 
setMethod('sparseness', signature(x='KINOMO'), 
	function(x){		
		# return the sparseness of the basis and coef matrix
		c(basis=sparseness(basis(x)), coef=sparseness(coef(x)))
	}
)



#' 
setGeneric('purity', function(x, y, ...) standardGeneric('purity') )
#' Computes the purity directly from the contingency table \code{x}
setMethod('purity', signature(x='table', y='missing'), 
	function(x, y){
		#for each cluster: compute maximum number of samples common to a class
		t <- apply(x, 1, max)
		# average and return the result
		sum(t) / sum(x)
	}
)
#' Computes the purity on the contingency table of \code{x} and \code{y}, that is  
#' coerced into a factor if necessary.
setMethod('purity', 'factor', 
	function(x, y, ...){
		
		# coerce `y` into a factor if necessary
		if( !is.factor(y) ) y <- as.factor(y)
		#compute the purity on the contingency table between clusters and true classes (clusters are in rows)
		purity(table(x, y), ...)
	}
)
#' Default method that should work for results of clustering algorithms, that have a 
#' suitable \code{predict} method that returns the cluster membership vector:
#' the purity is computed between \code{x} and \code{predict{y}}
setMethod('purity', 'ANY', 
	function(x, y, ...){
		# compute the purity for the samples clusters defined by the profiles
		purity(predict(x), y, ...)
	}
)


#' 
setGeneric('entropy', function(x, y, ...) standardGeneric('entropy') )
#' Computes the purity directly from the contingency table \code{x}.
#' 
#' This is the workhorse method that is eventually called by all other methods.
setMethod('entropy', signature(x='table', y='missing'), 
	function(x, y, ...){
		#for each cluster: compute the inner sum
		t <- apply(x, 1, function(n){ c.size <- sum(n); n %*% ifelse( n!=0, log2(n/c.size), 0)} )
		
		# weight and return the result
		- sum(t) / ( sum(x) * log2(ncol(x)) )
	}
)
#' Computes the purity on the contingency table of \code{x} and \code{y}, that is 
#' coerced into a factor if necessary.
setMethod('entropy', 'factor', 
	function(x, y, ...){
		# coerce `y` into a factor if necessary
		if( !is.factor(y) ) y <- as.factor(y)
		#copmute entropy on contingency table between clusters and true classes (clusters are in rows)
		entropy(table(x, y))
	}
)
#' Default method that should work for results of clustering algorithms, that have a 
#' suitable \code{predict} method that returns the cluster membership vector:
#' the purity is computed between \code{x} and \code{predict{y}}
setMethod('entropy', 'ANY', 
	function(x, y, ...){
		# compute the entropy for the samples clusters defined by the metagenes expression matrix
		entropy(predict(x), y)
	}
)




#setGeneric('KINOMOApply', function(object, ...) standardGeneric('KINOMOApply') )
KINOMOApply <- function(X, MARGIN, FUN, ..., simplify = TRUE, USE.NAMES = TRUE){
	if( MARGIN == 1L )
		apply(basis(X), 1L, FUN, ...)
	else if( MARGIN == 4L )
		apply(basis(X), 2L, FUN, ...)
	else if( MARGIN == 2L )
		apply(coef(X), 2L, FUN, ...)
	else if( MARGIN == 5L )
		apply(coef(X), 1L, FUN, ...)
	else if( MARGIN == 3L ){
		b <- basis(X)
		p <- coef(X)
		sapply(setNames(seq(nbasis(X), basisnames(X)))
				, function(i, ...) FUN(b[,i], p[i,], ...)
				, simplify = simplify, USE.NAMES = USE.NAMES)
	}else stop("invalid argument 'MARGIN' (expected values are: 1-basis rows, 2-coef columns, 3-(basis columns, coef rows), or 4-basis columns or 5-coef rows)")
}

###% Utility function to compute the dominant column for each row for a matrix.
.predict.KINOMO <- function(x, prob=FALSE){
	
	if( !is.matrix(x) ) stop('KINOMO:::.predict.KINOMO : only works on matrices')
	if( !prob ){
		#for each column return the (row) index of the maximum
		return( as.factor(apply(x, 1L, function(v) which.max(abs(v)))) )
	}
	else{
		#for each column return the (row) index of the maximum AND the associated probaility
		res <- apply(x, 1L,
				function(p){
					p <- abs(p)
					i <- which.max(p)
					c(i, p[i]/sum(p))
				}
		)
		# return the result as a list of two elements
		return( list(predict=as.factor(res[1,]), prob=res[2,]) )
	}
}




setGeneric('predict', package='stats')



setMethod('predict', 'KINOMO',
		function(object, what=c('columns', 'rows', 'samples', 'features'), prob=FALSE, dmatrix = FALSE){
			# determine which matrix to use for the prediction
			what <- match.arg(what)
			x <- if( what %in% c('features', 'rows') ) basis(object, all=FALSE) else t(coef(object, all=FALSE))
			
			# compute the indice of the dominant row for each column
            res <- .predict.KINOMO(x, prob)
            # attach dissimilarity matrix if requested
            if( dmatrix ){
                attr(res, 'dmatrix') <- 1 - cor(t(x))
            }
			return( res )
		}
)




#' 
setGeneric('basiscor', function(x, y, ...) standardGeneric('basiscor') )
#' Computes the correlations between the basis vectors of \code{x} and 
#' the columns of \code{y}. 
setMethod('basiscor', signature(x='KINOMO', y='matrix'),
	function(x, y, ...){
		cor(basis(x), y, ...)
	}
)
#' Computes the correlations between the columns of \code{x} 
#' and the the basis vectors of \code{y}.  
setMethod('basiscor', signature(x='matrix', y='KINOMO'),
	function(x, y, ...){
		cor(x, basis(y), ...)
	}
)
#' Computes the correlations between the basis vectors of \code{x} and \code{y}.
setMethod('basiscor', signature(x='KINOMO', y='KINOMO'),
		function(x, y, ...){
			basiscor(x, basis(y), ...)
		}
)
#' Computes the correlations between the basis vectors of \code{x}.
setMethod('basiscor', signature(x='KINOMO', y='missing'),
		function(x, y, ...){
			basiscor(x, x, ...)
		}
)



#' 
setGeneric('profcor', function(x, y, ...) standardGeneric('profcor') )
#' Computes the correlations between the basis profiles of \code{x} and 
#' the rows of \code{y}.
setMethod('profcor', signature(x='KINOMO', y='matrix'),
		function(x, y, ...){
			cor(t(coef(x)), t(y), ...)
		}
)
#' Computes the correlations between the rows of \code{x} and the basis 
#' profiles of \code{y}.
setMethod('profcor', signature(x='matrix', y='KINOMO'),
		function(x, y, ...){
			cor(t(x), t(coef(y)), ...)
		}
)
#' Computes the correlations between the basis profiles of \code{x} and \code{y}.
setMethod('profcor', signature(x='KINOMO', y='KINOMO'),
		function(x, y, ...){
			profcor(x, coef(y), ...)
		}
)
#' Computes the correlations between the basis profiles of \code{x}.
setMethod('profcor', signature(x='KINOMO', y='missing'),
		function(x, y, ...){
			profcor(x, x, ...)
		}
)



#' 
setGeneric('connectivity', function(object, ...) standardGeneric('connectivity') )


setMethod('connectivity', 'ANY', 
	function(object, ...){
		c <- predict(object, ...);
		outer(c, c, function(x,y) ifelse(x==y, 1,0));
	}
)



setMethod('connectivity', 'factor', 
	function(object, ...){
		outer(object, object, function(x,y) ifelse(x==y, 1,0));
	}
)
#' Equivalent to \code{connectivity(as.factor(x))}. 
setMethod('connectivity', 'numeric', 
	function(object, ...){
		connectivity(as.factor(object), ...)
	}
)


setMethod('connectivity', 'KINOMO', 
	function(object, no.attrib=FALSE){
		C <- callNextMethod(object=object, what='samples');
		if( !no.attrib ){
			class(C) <- c(class(C), 'KINOMO.consensus')
			attr(C, 'model') <- object
			attr(C, 'nrun') <- 1
			attr(C, 'nbasis') <- nbasis(object)
		}
		C
	}
)

# Unit test
unit.test(connectivity,{
	
	# build reference matrix
	n <- 10
	ref <- matrix(0, 2*n, 2*n)
	ref[1:n,1:n] <- 1
	ref[(n+1):(2*n),(n+1):(2*n)] <- 1
	
	checkIdentical(connectivity(gl(2, n)), ref, 'Factor')
	checkIdentical(connectivity(as.numeric(gl(2, n))), ref, 'Vector')
	# test with KINOMO model
	i <- gl(2, n)
	x <- KINOMOModel(H=matrix(c(rev(i), i), 2, byrow=TRUE))
	checkEquals(connectivity(x), ref, 'KINOMO model', check.attributes = FALSE)
	s <- sample.int(2*n)
	checkEquals(connectivity(x[,s]), ref[s,s], 'KINOMO model (shuffled)', check.attributes = FALSE)
	
})



setGeneric('rss', function(object, ...) standardGeneric('rss'))


setMethod('rss', 'matrix', 
	function(object, target){
		# make sure the target is provided
		if( missing(target) ) stop("KINOMO::rss - Argument 'target' is missing and required to compute the residual sum of squares.")
		
		# use the expression matrix if necessary
		if( inherits(target, 'ExpressionSet') ){
			# requires Biobase
			if( !require.quiet("Biobase") ) 
				stop("KINOMO::rss - The 'Biobase' package is required to extract expression data from 'ExpressionSet' objects [see ?'KINOMO-bioc']")
			
			target <- Biobase::exprs(target)
		}else if( is.data.frame(target) )
			target <- as.matrix(target)
		
		# return rss using the optimized C function
		.rss(object,target)
	}
)


setMethod('rss', 'ANY', 
	function(object, ...){
		rss(fitted(object), ...)
	}
)


unit.test(rss, {
	
	x <- rmatrix(20,10, max=50)
	y <- rmatrix(20,10, max=50)
	checkIdentical(rss(x, y), sum((x-y)^2), "Random matrices")
	
	y <- rKINOMO(3, x) # random compatible model
	r1 <- rss(y, x)
	checkIdentical(r, sum((x-fitted(y))^2), 'KINOMO model')
	checkIdentical(rss(y, Biobase::ExpressionSet(x)), sum((x-fitted(y))^2), 'KINOMO model (ExpressionSet)')
	y <- KINOMO(x, 3)
	r2 <- rss(y, x)
	checkIdentical(r2, sum((x-fitted(y))^2), 'Fitted KINOMO model')
	checkTrue(r2 < r1, 'Fitted KINOMO model has better RSS')
	y <- KINOMO(x, 3, 'lee')
	checkTrue(rss(y, x) < r2, "Fitted KINOMO model with 'lee' has better RSS than 'brunet'")
	
})


setGeneric('evar', function(object, ...) standardGeneric('evar'))


setMethod('evar', 'ANY', 
	function(object, target, ...){
		
		# make sure the target is provided
		if( missing(target) ) stop("KINOMO::evar - Argument 'target' is missing and required to compute the explained variance.")
		
		# use the expression matrix if necessary
		if( inherits(target, 'ExpressionSet') ){
			# requires Biobase
			if( !require.quiet("Biobase") ) 
				stop("KINOMO::evar - The 'Biobase' package is required to extract expression data from 'ExpressionSet' objects [see ?'KINOMO-bioc']")
			
			target <- Biobase::exprs(target)
		}
		
		t <- as.numeric(target)
		1 - rss(object, target, ...) / sum(t^2)
	}
)



setGeneric('deviance', package='stats')


#' 
setMethod('deviance', 'KINOMO', 
	function(object, y, method=c('', 'KL', 'euclidean'), ...){

	fun <- KINOMODistance(method)
	
	if( is.null(fun) ){
		warning('Undefined distance method: distance cannot be computed [returned NA]')
		return(as.numeric(NA))
	}
	
	# extract expression data from ExpressionSet objects
	if( is(y, 'ExpressionSet') )
		y <- Biobase::exprs(y)
	
	# apply the function and return the result
	fun(object, y, ...)

	}
)



KINOMODistance <- function(method=c('', 'KL', 'euclidean')){
	
		#message('compute distance')
		# determinate the distance measure to use
		if( is.null(method) ) return(NULL)
		
		if( is.character(method) ){
			errMeth <- try(method <- match.arg(method), silent=TRUE)
			# if the method is not predefined, try to find a function with the given name
			if( inherits(errMeth, 'try-error') ){			
				#TODO: this is not working with local functions
				if( is.character(method) ){
					errFun <- try(fun <- match.fun(method), silent=TRUE)
					if( inherits(errFun, 'try-error') ) 
						stop("Could not find distance measure '", method, "':\n\t- not a predefined measures -> ", errMeth,"\t- not a function -> ", errFun)
				}
				else fun <- method
				
				if( !is.function(fun) )
					stop('Invalid distance measure: should be a character string or a valid function definition')
			}
			else{
				# compute and return the distance measure		
				fun <- switch(method,
						euclidean = function(x, y, ...){
							# call optimized C function
							.rss(y, fitted(x))/2							
						},
						KL = function(x, y, ...){							
							# call optimized C function
							.KL(y, fitted(x))					
						}
				)
			}
		}
		else if( is.function(method) )
			fun <- method
		else
			stop('Invalid distance measure: should be a character string or a valid function definition')
	
		# return the distance function
		fun
	
	}
		




setGeneric('KINOMO.equal', function(x, y, ...) standardGeneric('KINOMO.equal') )
#' Compares two KINOMO models.
#' 
#' Arguments in \code{...} are used only when \code{identical=FALSE} and are 
#' passed to \code{all.equal}.
setMethod('KINOMO.equal', signature(x='KINOMO', y='KINOMO'), 
		function(x, y, identical=TRUE, ...){
			
			dots <- list(...)
			if( identical && length(dots) == 0 )
				identical(x, y)
			else
				all.equal(x, y, ...)
		}
)

# Match and Order Basis Components
#
# 
match.basis <- function(object, return.table=FALSE){
	
	# compute the contingency table
	#pcmap <- predict(object, 'cmap')
	# build the tree from consensus matrix
	h <- hclust(as.dist(1-consensus(object)), method='average')
	# extract membership from the tree
	cl <- cutree(h, k=nbasis(object))
	# change the class indexed to match the order of the consensus clusters 
	cl <- match(cl, unique(cl[h$order]))
	pcmap <- as.factor(cl)
	occ <- table(consensus=pcmap, fit=predict(object))		
	
	# add names if present
#	if( !is.null(basisnames(object)) ){
#		rownames(occ) <- colnames(occ) <- basisnames(object)
#	}
		
	# for each estimated component look for the maximum agreement		
	T.tmp <- occ				
	res <- rep(0, ncol(T.tmp))	
	for( i in 1:ncol(T.tmp) ){
		# get the row and column index of the maximum over the remaining entries 
		xm <- which.max(T.tmp)-1
		jm <- xm %/% nrow(T.tmp) + 1
		im <- xm - (jm-1) * nrow(T.tmp) + 1
		
		# assign the estimate row to the inferred reference column
		stopifnot( res[im]==0 )
		res[im] <- jm
		
		# erase the assigned estimate row
		T.tmp[im,] <- NA
		# erase the assigned reference column
		T.tmp[,jm] <- NA		
	}
	
	# return the mapping as an integer vector
	res <- as.integer(res)
	if( return.table )
		res <- list(match=res, table=occ)
	
	# return result
	res
}
