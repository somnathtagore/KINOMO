#' @include rmatrix.R
#' @include KINOMO-package.R
NULL

#' Utility Function in the KINOMO Package
#' 
#' @name utils-KINOMO
#' @rdname utils
NULL

#' Internal verbosity option
#' @param val logical that sets the verbosity level.
#' @return the old verbose level   
#' @keywords internal
lverbose <- local({
			.val <- 0
			function(val){
				if( missing(val) ) return(.val)
				oval <- .val
				.val <<- val
				invisible(oval)
			}
		})
vmessage <- function(..., appendLF=TRUE) if( lverbose() ) cat(..., if(appendLF) "\n", sep='')


KINOMO_stop <- function(name, ...){
	stop("KINOMO::", name , ' - ', ..., call.=FALSE)
}
KINOMO_warning <- function(name, ...){
	warning("KINOMO::", name , ' - ', ..., call.=FALSE)
}

# or-NULL operator
'%||%' <- function(x, y) if( !is.null(x) ) x else y 

# cat object or class for nice cat/message
quick_str <- function(x) if( is.atomic(x) ) x else class(x)[1]

# remove all attributes from an object
rmAttributes <- function(x){
	attributes(x) <- NULL
	x
}


str_args <- function(x, exdent=10L){
	s <- capture.output(print(args(x)))
	paste(str_trim(s[-length(s)]), collapse=str_c('\n', paste(rep(' ', exdent), collapse='')))
}


txtProgressBar <- function (min = 0, max = 1, initial = 0, char = "=", width = NA, 
    	title= if( style == 3 ) ' ', label, style = 1, file = ""
		, shared = NULL)
{
    if (!identical(file, "") && !(inherits(file, "connection") && 
        		isOpen(file))) 
        stop("'file' must be \"\" or an open connection object")
    if (!style %in% 1L:3L) 
        style <- 1
    .val <- initial
    .killed <- FALSE
    .nb <- 0L
    .pc <- -1L
    nw <- nchar(char, "w")
    if (is.na(width)) {
        width <- getOption("width")
        if (style == 3L) 
            width <- width - 10L
        width <- trunc(width/nw)
    }
    if (max <= min) 
        stop("must have max > min")
	
	# setup shared directory
	.shared <- NULL
	if( isTRUE(shared) ) shared <- tempdir()
	if( is.character(shared) ){
		.shared <- tempfile('progressbar_', tmpdir=shared[1L])
		dir.create(.shared)
	}
	#
	
	getval <- function(value){
		if( value >= max || value <= min 
			|| is.null(.shared) || !file.exists(.shared) ){
			value
		}else{
			cat('', file=file.path(.shared, paste('_', value, sep='')))
			length(list.files(.shared))
		}
	}
	
    up1 <- function(value) {
		if (!is.finite(value) || value < min || value > max) 
            return()
		
		# get actual value 
		value <- getval(value)
		
        .val <<- value
        nb <- round(width * (value - min)/(max - min))
        if (.nb < nb) {
            cat(paste(rep.int(char, nb - .nb), collapse = ""), 
                	file = file)
            flush.console()
        }
        else if (.nb > nb) {
            cat("\r", title, paste(rep.int(" ", .nb * nw), collapse = ""), 
                	"\r", title, paste(rep.int(char, nb), collapse = ""), 
                	sep = "", file = file)
            flush.console()
        }
        .nb <<- nb
    }
    up2 <- function(value) {
        if (!is.finite(value) || value < min || value > max) 
            return()
		
		# get actual value 
		value <- getval(value)
		
        .val <<- value
        nb <- round(width * (value - min)/(max - min))
        if (.nb <= nb) {
            cat("\r", title, paste(rep.int(char, nb), collapse = ""), 
                	sep = "", file = file)
            flush.console()
        }
        else {
            cat("\r", title, paste(rep.int(" ", .nb * nw), collapse = ""), 
                	"\r", paste(rep.int(char, nb), collapse = ""), 
                	sep = "", file = file)
            flush.console()
        }
        .nb <<- nb
    }
    up3 <- function(value) {
        if (!is.finite(value) || value < min || value > max) 
            return()
		
		# get actual value 
		value <- getval(value)
		
        .val <<- value
        nb <- round(width * (value - min)/(max - min))
        pc <- round(100 * (value - min)/(max - min))
        if (nb == .nb && pc == .pc) 
            return()
        cat(paste(c("\r",title," |", rep.int(" ", nw * width + 6)), collapse = ""), 
            	file = file)
        cat(paste(c("\r",title," |", rep.int(char, nb), rep.int(" ", 
            							nw * (width - nb)), sprintf("| %3d%%", pc)), collapse = ""), 
            	file = file)
        flush.console()
        .nb <<- nb
        .pc <<- pc
    }
    getVal <- function() .val
    kill <- function(cleanup=TRUE) if (!.killed) {
        	cat("\n", file = file)
        	flush.console()
			.killed <<- TRUE
			
			# do some cleanup
			if( cleanup ){
				# delete shared directory
				if( !is.null(.shared) && file.exists(.shared) ) 
					unlink(.shared, recursive=TRUE)
				#
			}
			invisible(TRUE)
    	}
    up <- switch(style, up1, up2, up3)
    up(initial)
    structure(list(getVal = getVal, up = up, kill = kill), class = "txtProgressBar")
}

###% apply a function to each entry in a matrix
matapply <- function(x, FUN, ...){
	res <- sapply(x, FUN, ...)
	matrix(res, nrow(x))
}

###% try to convert a character string into a numeric
toNumeric <- function(x){
	suppressWarnings( as.numeric(x) )
}

###% Tells one is running in Sweave
isSweave <- function() !is.null(sweaveLabel()) 
	
sweaveLabel <- function(){
	if ((n.parents <- length(sys.parents())) >= 3) {
		for (i in seq_len(n.parents) - 1) {
			if ("chunkopts" %in% ls(envir = sys.frame(i))) {
				chunkopts = get("chunkopts", envir = sys.frame(i))
				if (all(c("prefix.string", "label") %in% names(chunkopts))) {
					img.name = paste(chunkopts$prefix.string, chunkopts$label, 
							sep = "-")
					return(img.name)
					break
				}
			}
		}
	}
}

sweaveFile <- function(){
	label <- sweaveLabel()
	if( !is.null(label) )
		paste(label, '.pdf', sep='')
}

fixSweaveFigure <- function(filename){
	if( missing(filename) ){
		filename <- sweaveLabel()
		if( is.null(filename) ) return()
		filename <- paste(filename, '.pdf', sep='')
	}
	filepath <- normalizePath(filename)
	tf <- tempfile()
	system(paste("pdftk", filepath, "cat 2-end output", tf, "; mv -f", tf, filepath))
}

###% 'more' functionality to read data progressively
more <- function(x, step.size=10, width=20, header=FALSE, pattern=NULL){
	
	if( !(is.matrix(x) || is.data.frame(x) || is.vector(x) || is.list(x)) )
		stop("KINOMO::more - invalid argument 'x': only 'matrix', 'data.frame', 'vector' and 'list' objects are handled.")
	
	one.dim <- is.null(dim(x))
	single.char <- FALSE
	n <-
		if( is.character(x) && length(x) == 1 ){			
			cat("<character string:", nchar(x), ">\n")
			single.char <- TRUE
			nchar(x)
		}
		else if( one.dim ){
			cat("<", class(x),":", length(x), ">\n")
			
			# limit to matching terms if necessary
			if( !is.null(pattern) )
				x[grep(pattern, x)]
			
			length(x)
		}else{
			cat("<", class(x),":", nrow(x), "x", ncol(x), ">\n")
			head.init <- colnames(x)
			head.on <- TRUE
			
			# limit to matching terms if necessary
			if( !is.null(pattern) ){
				idx <- apply(x, 2, grep, pattern=pattern)
				print(idx)
				idx <- unique(if( is.list(idx) ) unlist(idx) else as.vector(idx))
				x <- x[idx,, drop=FALSE]
			}
			
			nrow(x)
		}	
		
	i <- 0
	while( i < n ){
		# reduce 'step.size' if necessary
		step.size <- min(step.size, n-i)
		
		what2show <- if( single.char )
			substr(x, i+1, i+step.size)
		else if( one.dim )			
			if( !is.na(width) ) sapply(x[seq(i+1, i+step.size)], function(s) substr(s, 1, width) ) else x[seq(i+1, i+step.size)]
		else{
			w <- x[seq(i+1, i+step.size), , drop=FALSE]
			if( !is.na(width) ){ 
				w <- apply(w, 2, 
					function(s){
						ns <- toNumeric(s)
						if( !is.na(ns[1]) ) # keep numerical value as is
							ns
						else # limit output if required
							substr(s, 1, width)
						
					}) 
				rownames(w) <- rownames(x)[seq(i+1, i+step.size)]
			} 
				
			
			# remove header if not required
			if( !header && head.on ){
				colnames(x) <- sapply(colnames(x), function(c) paste(rep(' ', nchar(c)), collapse=''))
				head.on <- FALSE
			}
			
			# return the content
			w
		}
		
		cat( show(what2show) )
		i <- i + step.size
		
		# early break if necessary
		if( i >= n )
			break
		# ask user what to to next
		ans <- scan(what='character', quiet=TRUE, n=1, multi.line=FALSE)
		
		# process user command if any (otherwise carry on)
		if( length(ans) > 0 ){		
			if( !is.na(s <- toNumeric(ans)) ) # change step size
				step.size <- s
			else if( !header && ans %in% c('h', 'head') ){
				colnames(x) <- head.init
				head.on <- TRUE
			}
			else if( ans %in% c('q', 'quit') ) # quit
				break
		}
	}
	invisible()
} 


#' 
randomize <- function(x, ...){
	
	if( is(x, 'ExpressionSet') ) x <- Biobase::exprs(x)
		
	# resample the columns
	apply(x, 2, function(c, ...) sample(c, size=length(c), ...), ...)
	
}

###% Returns the rank-k truncated SVD approximation of x
tsvd <- function(x, r, ...){
	stopifnot( r > 0 && r <= min(dim(x)))
	s <- svd(x, nu=r, nv=r, ...)
	s$d <- s$d[1:r]
	
	# return results
	s
}

###% Subset a list leaving only the arguments from a given function 
.extract.args <- function(x, fun, ...){
	
	fdef <- if( is.character(fun) )	getFunction(fun, ...)
			else if( is.function(fun) ) fun
			else stop("invalid argument 'fun': expected function name or definition")
	
	if( length(x) == 0 ) return(x)	
	x.ind <- charmatch(if( is.list(x) ) names(x) else x, args <- formalArgs(fdef))
	x[!is.na(x.ind)]
}

###% Returns the version of the package
KINOMOInfo <- function(command){	
	pkg <- 'KINOMO'
	curWarn <- getOption("warn")
	on.exit(options(warn = curWarn), add = TRUE)
	options(warn = -1)
	desc <- packageDescription(pkg, fields="Version")
	if (is.na(desc)) 
		stop(paste("Package", pkg, "not found"))
	desc
}

###% Returns TRUE if running under Mac OS X + GUI
is.Mac <- function(check.gui=FALSE){
	is.mac <- (length(grep("darwin", R.version$platform)) > 0)
	# return TRUE is running on Mac (adn optionally through GUI)
	is.mac && (!check.gui || .Platform$GUI == 'AQUA')
}

###% Hash a function body (using digest)
#' @import digest
hash_function <- function(f){
	b <- body(f)
	attributes(b) <- NULL
	fdef <- paste(c(capture.output(args(f))[1], capture.output(print(b))), collapse="\n")
	# print(fdef)
	digest(b)
}


###% compare function with copy and with no copy
cmp.cp <- function(...){
	res <- KINOMO(..., copy=F)
	resc <- KINOMO(..., copy=T)
	cat("identical: ", identical(fit(res), fit(resc))
			, " - all.equal: ", all.equal(fit(res), fit(resc))
			, " - diff: ", all.equal(fit(res), fit(resc), tol=0)
			, "\n"
	)
	invisible(res)
} 

# return the internal pointer address 
C.ptr <- function(x, rec=FALSE)
{	
	attribs <- attributes(x)
	if( !rec || is.null(attribs) )
		.Call("ptr_address", x, PACKAGE='KINOMO')
	else
		c( C.ptr(x), sapply(attribs, C.ptr, rec=TRUE))
	
}

is.same <- function(x, y){
	C.ptr(x) == C.ptr(y)
}

is.eset <- function(x) is(x, 'ExpressionSet')

# clone an object
clone <- function(x){
	.Call('clone_object', x, PACKAGE='KINOMO')
}

# deep-clone an object
clone2 <- function(x){
	if( is.environment(x) ){
		y <- Biobase::copyEnv(x)
		eapply(ls(x, all.names=TRUE), 
			function(n){
				if( is.environment(x[[n]]) ){
					y[[n]] <<- clone(x[[n]])
					if( identical(parent.env(x[[n]]), x) )
						parent.env(y[[n]]) <<- y
				}
		})
	}else{
		y <- .Call('clone_object', x, PACKAGE='KINOMO')		
		if( isS4(x) ){ ## deep copy R object
			lapply(slotNames(class(y)), 
				function(n){					
					slot(y, n) <<- clone(slot(x, n)) 
			})
		}else if( is.list(x) ){ ## copy list or vector
			sapply(seq_along(x), 
				function(i){					
					y[[i]] <<- clone(x[[i]])					
			})
		}
	}
	
	y
}

#compute RSS with C function
.rss <- function(x, y)
{	
	.Call("Euclidean_rss", x, y, PACKAGE='KINOMO')
}

#compute KL divergence with C function
.KL <- function(x, y)
{	
	.Call("KL_divergence", x, y, PACKAGE='KINOMO')
}


pmax.inplace <- function(x, lim, skip=NULL){
	
	.Call('ptr_pmax', x, lim, as.integer(skip), PACKAGE='KINOMO')
	
}

# colMin
colMin <- function(x){
	.Call('colMin', x, PACKAGE='KINOMO')
}

# colMax
colMax <- function(x){
	.Call('colMax', x, PACKAGE='KINOMO')
}


neq.constraints.inplace <- function(x, constraints, ratio=NULL, value=NULL, copy=FALSE){
	
	# if requested: clone data as neq.constrains.inplace modify the input data in place
	if( copy )
		x <- clone(x)
	
	.Call('ptr_neq_constraints', x, constraints, ratio, value, PACKAGE='KINOMO')	
}

# Test if an external pointer is nil
# Taken from package bigmemory 
ptr_isnil <- function (address) 
{
	if (class(address) != "externalptr") 
		stop("address is not an externalptr.")
	.Call("ptr_isnil", address, PACKAGE='KINOMO')	
}


###% Draw the palette of colors
###% 
###% Taken from the examples of colorspace::rainbow_hcl
###% 
pal <- function(col, h=1, border = "light gray")
{
	n <- length(col)	
	plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, h), axes = FALSE, xlab = "", ylab = "")
	rect(0:(n-1)/n, 0, 1:n/n, h, col = col, border = border)
}
###% Draw the Palette of Colors as a Wheel
###% 
###% Taken from the examples of colorspace::rainbow_hcl
###% 
wheel <- function(col, radius = 1, ...)
	pie(rep(1, length(col)), col = col, radius = radius, ...)

# Define a S4 class to handle function slots given as either a function definition 
# or a character string that gives the function's name. 
setClassUnion('.functionSlot', c('character', 'function'))

# Define a S4 class to handle function slots given as either a function definition 
# or a character string that gives the function's name or NULL.
setClassUnion('.functionSlotNULL', c('character', 'function', 'NULL'))
.validFunctionSlot <- function(slot, allow.empty=FALSE, allow.null=TRUE){
	if( is.null(slot) ){
		if( !allow.null ) return('NULL value is not allowed')
		return(TRUE)
	}
	if( is.character(slot) ){
		if( !allow.empty && slot == '' ) return('character string cannot be empty')
		if( length(slot) != 1 ) return(paste('character string must be a single value [length =', length(slot), ']', sep=''))
	}			
	
	return(TRUE)
}



# Extracted from the psychometric package (0.1.0)
# Copyright Thomas D. Fletcher
# Under Gnu GPL2
CI.Rsqlm <- function (obj, level = 0.95) 
{
	l <- level
	rsq <- summary(obj)$r.squared
	k <- summary(obj)$df[1] - 1
	n <- obj$df + k + 1
	mat <- CI.Rsq(rsq, n, k, level = l)
	return(mat)
}
# Extracted from the psychometric package (0.1.0)
# Copyright Thomas D. Fletcher
# Under Gnu GPL2
CI.Rsq <- function (rsq, n, k, level = 0.95) 
{
	noma <- 1 - level
	sersq <- sqrt((4 * rsq * (1 - rsq)^2 * (n - k - 1)^2)/((n^2 - 
							1) * (n + 3)))
	zs <- -qnorm(noma/2)
	mez <- zs * sersq
	lcl <- rsq - mez
	ucl <- rsq + mez
	mat <- data.frame(Rsq = rsq, SErsq = sersq, LCL = lcl, UCL = ucl)
	return(mat)
}

str_dim <- function(x, dims=dim(x)){
    if( !is.null(dims) ) paste0(dims, collapse = ' x ')
    else length(x)
}

# Internal override stringr function str_match
# 
# This is to get the previous behaviour on optional groups, because 
# in stringr >= 1.0.0 absent optional groups get an NA value instead 
# of an empty string, which in turn breaks some downstream processing.
str_match <- function(...){
    
    res <- stringr::str_match(...)
    # replace NAs by "" for globally matched strings
    if( length(w <- which(!is.na(res[, 1L]))) ){
        res[w, ][is.na(res[w, ])] <- ""
    }
    res
}

