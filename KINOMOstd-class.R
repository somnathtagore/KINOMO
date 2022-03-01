

#' @include KINOMO-class.R
NULL





#'  
setClass('KINOMOstd' 
		, representation(
			W = 'matrix' # basis matrix
			, H = 'matrix' # mixture coefficients matrix
			, bterms = 'data.frame' # fixed basis terms: nrow(bterms) = nrow(x)
			, ibterms = 'integer' # index of the fixed basis terms
			, cterms = 'data.frame' # fixed coef terms: ncol(cterms) = ncol(x)
			, icterms = 'integer' # index of the fixed coefficient terms
		)
		
		, prototype = prototype(
				W = matrix(as.numeric(NA), 0, 0),
				H = matrix(as.numeric(NA), 0, 0)
		)
		
		, validity = function(object){
			
			# dimension compatibility: W and H must be compatible for matrix multiplication
			if( ncol(object@W) != nrow(object@H) ){
				return(paste('Dimensions of W and H are not compatible [ncol(W)=', ncol(object@W) , '!= nrow(H)=', nrow(object@H), ']'))
			}
			# give a warning if the dimensions look strange: rank greater than the number of samples
			if( !is.empty.KINOMO(object) && ncol(object@H) && ncol(object@W) > ncol(object@H) ){
				warning(paste('Dimensions of W and H look strange [ncol(W)=', ncol(object@W) , '> ncol(H)=', ncol(object@H), ']'))
			}
			
			# everything went fine: return TRUE
			return(TRUE)
		}
		, contains = 'KINOMO'
)



#' 
setMethod('.basis', 'KINOMOstd',
	function(object){ 
		object@W
	}
)
#' Set the basis matrix in standard KINOMO models 
#' 
#' This function sets slot \code{W} of \code{object}.
setReplaceMethod('.basis', signature(object='KINOMOstd', value='matrix'), 
	function(object, value){ 
		object@W <- value		
		object
	} 
)

#' Get the mixture coefficient matrix in standard KINOMO models 
#' 
#' This function returns slot \code{H} of \code{object}.
setMethod('.coef', 'KINOMOstd',
	function(object){
		object@H
	}
)
#' Set the mixture coefficient matrix in standard KINOMO models 
#' 
#' This function sets slot \code{H} of \code{object}.
setReplaceMethod('.coef', signature(object='KINOMOstd', value='matrix'), 
	function(object, value){ 
		object@H <- value			
		object
	}
)



#' 
setMethod('fitted', signature(object='KINOMOstd'), 
	function(object, W, H, ...){
		if( missing(W) ) W <- object@W
		if( missing(H) ) H <- object@H
		return(W %*% H)
	}
)

