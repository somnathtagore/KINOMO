

#' @include KINOMO-class.R
NULL




#' 
setGeneric('featureScore', function(object, ...) standardGeneric('featureScore') )
#' Computes feature scores on a given matrix, that contains the basis component in columns. 
setMethod('featureScore', 'matrix', 
	function(object, method=c('kim', 'max')){
		
		method <- match.arg(method)
		score <- switch(method,
				
				kim = {
					#for each row compute the score
					s <- apply(object, 1, function(g){
								g <- abs(g)
								p_i <- g/sum(g)
								crossprod(p_i, log2(p_i))
							})
					# scale, translate and return the result
					1 + s / log2(ncol(object))		
				}
				, max = {
					apply(object, 1L, function(x) max(abs(x)))
				}
		)
		
		# return the computed score
		return(score)
	}
)
#' Computes feature scores on the basis matrix of an KINOMO model.
setMethod('featureScore', 'KINOMO', 
	function(object, ...){
		featureScore(basis(object), ...)
	}
)



#' 
setGeneric('extractFeatures', function(object, ...) standardGeneric('extractFeatures') )

# internal functio to trick extractFeatures when format='subset'
.extractFeaturesObject <- local({
	.object <- NULL
	function(object){
		# first call resets .object
		if( missing(object) ){
			res <- .object
			.object <<- NULL
			res
		}else # set .object for next call
			.object <<- object
	}
})



#' 
setMethod('extractFeatures', 'matrix', 
	function(object, method=c('kim', 'max')
			, format=c('list', 'combine', 'subset'), nodups=TRUE){
		
		res <-
				if( is.numeric(method) ){
					# repeat single values
					if( length(method) == 1L ) method <- rep(method, ncol(object))
					
					# float means percentage, integer means count
					# => convert into an integer if values > 1						
					if( all(method > 1L) ) method <- as.integer(method)
					
					if( is.integer(method) ){ # extract top features
						
						# only keep the specified number of feature for each column
						mapply(function(i, l)	head(order(object[,i], decreasing=TRUE), l)
								, seq(ncol(object)), method, SIMPLIFY=FALSE)
						
					}else{ # extract features with contribution > threshold
						
						# compute relative contribution
						so <- sweep(object, 1L, rowSums(object), '/')
						# only keep features above threshold for each column
						mapply(function(i, l)	which(so[,i] >= l)
								, seq(ncol(object)), method, SIMPLIFY=FALSE)
						
					}
				}else{
					method <- match.arg(method)
					switch(method,
							kim = { # KIM & PARK method
								
								# first score the genes
								s <- featureScore(object, method='kim')
								
								# filter for the genes whose score is greater than \mu + 3 \sigma
								th <- median(s) + 3 * mad(s)
								sel <- s >= th
								#print( s[sel] )
								#print(sum(sel))
								
								# build a matrix with:
								#-> row#1=max column index, row#2=max value in row, row#3=row index
								temp <- 0;
								g.mx <- apply(object, 1L, 
										function(x){
											temp <<- temp +1
											i <- which.max(abs(x));
											#i <- sample(c(1,2), 1)
											c(i, x[i], temp)
										}
								)
								
								# test the second criteria
								med <- median(abs(object))
								sel2 <- g.mx[2,] >= med
								#print(sum(sel2))
								
								# subset the indices
								g.mx <- g.mx[, sel & sel2, drop=FALSE]				
								# order by decreasing score
								g.mx <- g.mx[,order(s[sel & sel2], decreasing=TRUE)]
								
								# return the indexes of the features that fullfil both criteria
								cl <- factor(g.mx[1,], levels=seq(ncol(object))) 
								res <- split(g.mx[3,], cl)
								
								# add the threshold used
								attr(res, 'threshold') <- th
								
								# return result
								res
								
							},
							max = { # MAX method from bioKINOMO
								
								# determine the specific genes for each basis vector
								res <- lapply(1:ncol(object), 
										function(i){
											mat <- object
											vect <- mat[,i]
											#order by decreasing contribution to factor i
											index.sort <- order(vect, decreasing=TRUE)		
											
											for( k in seq_along(index.sort) )
											{
												index <- index.sort[k]
												#if the feature contributes more to any other factor then return the features above it
												if( any(mat[index,-i] >= vect[index]) )
												{
													if( k == 1 ) return(as.integer(NA))
													else return( index.sort[1:(k-1)] )
												}
											}
											
											# all features meet the criteria
											seq_along(vect)
										}
								)
								
								# return res
								res
							}
					)
				}
		
		#Note: make sure there is an element per basis (possibly NA)
		res <- lapply(res, function(ind){ if(length(ind)==0) ind<-NA; as.integer(ind)} )
		
		# add names if possible
		if( !is.null(rownames(object)) ){
			noNA <- sapply(res, is_NA)
			res[noNA] <- lapply(res[noNA], function(x){
						setNames(x, rownames(object)[x])
					})
		}
		
		# apply the desired output format
		format <- match.arg(format)
		res <- switch(format
				#combine: return all the indices in a single vector
				, combine = { 
					# ensure that there is no names: for unlist no to mess up feature names
					names(res) <- NULL
					ind <- na.omit(unlist(res))
					if( nodups ) unique(ind) 
					else ind
				} 
				#subset: return the object subset with the selected indices
				, subset = {
					ind <- na.omit(unique(unlist(res)))
					sobject <- .extractFeaturesObject()
					{if( is.null(sobject) ) object else sobject}[ind, , drop=FALSE]
				}
				#else: leave as a list
				,{
					# add component names if any
					names(res) <- colnames(object)
					res
				}
		)
		
		# add attribute method to track the method used
		attr(res, 'method') <- method
		# return result
		return( res )
	}
)
#' Select basis-specific features from an KINOMO model, by applying the method 
#' \code{extractFeatures,matrix} to its basis matrix.
#' 
#'  
setMethod('extractFeatures', 'KINOMO', 
	function(object, ...){
		# extract features from the basis matrix, but subset the KINOMO model itself
		.extractFeaturesObject(object)
		extractFeatures(basis(object), ...)
	}
)

unit.test(extractFeatures, {
			
	.check <- function(x){
		msg <- function(...) paste(class(x), ':', ...)
		checkTrue( is.list(extractFeatures(x)), msg("default returns list"))
		checkTrue( is.list(extractFeatures(x, format='list')), msg("format='list' returns list"))
		checkTrue( is.integer(extractFeatures(x, format='combine')), msg("format='combine' returns an integer vector"))
		checkTrue( is(extractFeatures(x, format='subset'), class(x)), msg("format='subset' returns same class as object"))
	}
	
	.check(rmatrix(50, 5))
	.check(rKINOMO(3, 50, 5))
	
})
