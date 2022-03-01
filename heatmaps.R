

#' @include KINOMO-class.R
#' @include aheatmap.R
NULL


setGeneric('metaHeatmap', function(object, ...) standardGeneric('metaHeatmap') )

setMethod('metaHeatmap', signature(object='matrix'),
		function(object, ...){
			local <- function(object, type=c('plain', 'consensus'), class
				, unit.scaling=c('none', 'row', 'column'), palette="YlOrRd"
				, rev.palette=FALSE, show.prediction=TRUE, ...){
			
			.Defunct('metaHeatmap', 'KINOMO', "The S4 method 'metaHeatmap,matrix' is defunct, use 'aheatmap' instead.")
			




			}
			local(object, ...)
		}
)
#' Deprecated method that is substituted by \code{\link{coefmap}} and \code{\link{basismap}}.
setMethod('metaHeatmap', signature(object='KINOMO'),
		function(object, ...){
			local <- function(object, what=c('samples', 'features'), filter=FALSE, ...){
			
				what <- match.arg(what)
				if( what == 'samples' ){
					# send deprecated warning
					.Defunct('coefmap', 'KINOMO', "Direct use of the S4-Method 'metaHeatmap' for 'KINOMO' objects is defunct, use 'coefmap' instead.")
					
					# call the new function 'coefmap'
					return( coefmap(object, ...) )			
					
				}else if( what == 'features' ){
					# send deprecated warning
					.Defunct('basismap', 'KINOMO', "Direct use of the S4-Method 'metaHeatmap' for 'KINOMO' objects is defunct, use 'basismap' instead.")
					
					# call the new function 'basismap'
					return( basismap(object, subsetRow=filter, ...) )
					
				}
			}
			local(object, ...)
	}
)

# match an annotation track against list of supported tracks
match_named_track <- function(annotation, tracks, msg, optional=FALSE){
	
	idx <- 
	if( is.character(annotation) ){
		i <- match(annotation, tracks, nomatch=if(optional) 0L else NA )
		if( any(!is.na(i)) ){
			if( !optional && any(is.na(i)) ){
				stop(msg, "invalid track(s) [", str_out(annotation[is.na(i)])
						, "]: should be one of ", str_out(tracks))
			}
		}
		i
	}else if( is.list(annotation) ){ 
		sapply(annotation, function(x){
					if( isString(x) ) match(x, tracks, nomatch=if(optional) 0L else NA )
					else NA
				})
	}
	
	if( is.null(idx) ) return()
	ok <- !is.na(idx)
	# result
	# remaining annotations
	ann <- annotation[!ok]
	if( length(ann) == 0L ) ann <- NULL
	# track annotations
	tr <- unlist(annotation[which(ok)])
	idx <- idx[which(ok)] 
	if( is.null(names(annotation)) ) names(tr) <- tr
	else{
		mn <- names(tr) == ''
		names(tr)[mn] <- tr[mn]
	}
	others <- tr[idx==0L]
	#
#	list(ann=ann, tracks=tr[idx>0L], others=if(length(others)) others else NULL)
	list(ann=as.list(ann), tracks=tr)
}



#'  
NULL
 

setGeneric('basismap', function(object, ...) standardGeneric('basismap') )

setMethod('basismap', signature(object='KINOMO'),
	function(object, color = 'YlOrRd:50'
			, scale = 'r1' 
			, Rowv=TRUE, Colv=NA, subsetRow=FALSE
			, annRow=NA, annCol=NA, tracks = 'basis'
			, main="Basis components", info = FALSE
			, ...){
		
		# resolve subsetRow if its a single value
		if( is.atomic(subsetRow) && length(subsetRow) == 1 ){
			subsetRow <- 
				if( isFALSE(subsetRow) )
					NULL
				else if( isTRUE(subsetRow) ) # use Kim and Park scoring scheme for filtering 			
					extractFeatures(object, format='combine')
				else if( is.character(subsetRow) || is.numeric(subsetRow) ) # use subsetRow as a filtering method
					extractFeatures(object, method=subsetRow, format='combine')
				else stop("KINOMO::basismap - invalid single value for argument 'subsetRow' [logical, numeric or character expected]")
		}
		
		# extract the basis vector matrix
		x <- basis(object)
		
		# add side information if requested
		info <- if( isTRUE(info) && isKINOMOfit(object) ) 
					paste("Method:", algorithm(object))
				else if( isFALSE(info) ) NULL
				else info
		
		# process annotation tracks
		ptracks <- process_tracks(x, tracks, annRow, annCol)
		annRow <- ptracks$row 
		annCol <- ptracks$col
		# set special annotation handler
		specialAnnotation(1L, 'basis', function() predict(object, what='features'))
		specialAnnotation(2L, 'basis', function() as.factor(1:nbasis(object)))
		#
			
		# call aheatmap on matrix
		aheatmap(x, color = color, ...
				, scale = scale, Rowv=Rowv, Colv = Colv, subsetRow = subsetRow
				, annRow = annRow, annCol = annCol
				, main = main, info = info)	
	}
)

# check if an object contains some value
anyValue <- function(x){
	length(x) > 0L && !is_NA(x) 
}

grep_track <- function(x){
	list(
		both = grepl("^[^:].*[^:]$", x) | grepl("^:.*:$", x)
		, row = grepl("^:.*[^:]$", x)
		, col = grepl("^[^:].*:$", x)
	)
}

# process extra annotation tracks
process_tracks <- function(data, tracks, annRow=NA, annCol=NA){
	
	if( anyValue(tracks) ){
		
		# extract choices from caller function
		formal.args <- formals(sys.function(sys.parent()))
		choices <- eval(formal.args[[deparse(substitute(tracks))]])
		if( isTRUE(tracks) ) tracks <- choices
		else{
			if( !is.character(tracks) ) 
				stop("Special annotation tracks must be specified either as NA, TRUE or a character vector [", class(tracks), "].")
			
			# check validity
			pattern <- "^(:)?([^:]*)(:)?$"
			basech <- str_match(choices, pattern)
			basetr <- str_match(tracks, pattern)
			tr <- basetr[, 3L]
	#		print(basetr)
	#		print(basech)
			# extend base track name
			i <- charmatch(tr, basech[,3L])
			tr[!is.na(i)] <- basech[i[!is.na(i)],3L]
			tracks_long <- str_c(basetr[,2L], tr, basetr[,4L])
			# extend choices
			tty_choice <- grep_track(choices)
			if( any(tty_choice$both) )
				choices <- c(choices, str_c(':', choices[tty_choice$both]), str_c(choices[tty_choice$both], ':'))
			# look for exact match
			itr <- charmatch(tracks_long, choices)
			if( length(err <- which(is.na(itr))) ){
				stop("Invalid special annotation track name [", str_out(tracks[err], Inf)
					,"]. Should partially match one of ", str_out(choices, Inf), '.')
			}
			tracks[!is.na(itr)] <- choices[itr]
		}
#		print(tracks)
	}
	#
	tty <- grep_track(tracks)
	# create result object
	build <- function(x, ann, data, margin){
		t <- 
		if( anyValue(x) ) as.list(setNames(str_c(':', sub("(^:)|(:$)","",x)), names(x)))
		else NA
		# build annotations
		atrack(ann, t, .DATA=amargin(data,margin))
	}
	
	res <- list()
	res$row <- build(tracks[tty$both | tty$row], annRow, data, 1L)
	res$col <- build(tracks[tty$both | tty$col], annCol, data, 2L)
	#str(res)
	res
}


setGeneric('coefmap', function(object, ...) standardGeneric('coefmap') )


setMethod('coefmap', signature(object='KINOMO'),
		function(object, color = 'YlOrRd:50'
				, scale = 'c1'
				, Rowv = NA, Colv = TRUE
				, annRow = NA, annCol = NA, tracks='basis'
				, main="Mixture coefficients", info = FALSE
				, ...){
						
			# use the mixture coefficient matrix
			x <- coef(object)
			
			# add side information if requested
			info <- if( isTRUE(info) && isKINOMOfit(object) ) 
						paste("Method: ", algorithm(object))
					else if( isFALSE(info) ) NULL
					else info
			
			# process annotation tracks
			ptracks <- process_tracks(x, tracks, annRow, annCol)
			annRow <- ptracks$row
			annCol <- ptracks$col
			# set special annotation handler
			specialAnnotation(1L, 'basis', function() as.factor(1:nbasis(object)))
			specialAnnotation(2L, 'basis', function() predict(object))
			#
			
			## process ordering
			if( isString(Colv) ){
				if( Colv == 'basis' ) Colv <- 'samples'
				if( Colv == 'samples' )
					Colv <- order(as.numeric(predict(object, Colv)))
			}
			##
			
			# call aheatmap on matrix
			aheatmap(x, ..., color = color
					, scale = scale, Rowv = Rowv, Colv=Colv
					, annRow=annRow, annCol = annCol
					, main=main, info = info)
		}
)

