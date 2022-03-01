
###############################################################################

#' @include KINOMO-class.R
#' @include transforms.R
NULL



NULL



if(!requireNamespace("Biobase")) BiocManager::install("Biobase")

.onLoad.KINOMO.bioc <- function(){
	
if( pkgmaker::require.quiet('Biobase') ){

	# load Biobase package
	requireNamespace('Biobase')
	#library(Biobase)

	#' Performs KINOMO on an ExpressionSet: the target matrix is the expression matrix \code{exprs(x)}.
	#' @rdname bioc
	setMethod('KINOMO', signature(x='ExpressionSet', rank='ANY', method='ANY'), 
		function(x, rank, method, ...)
		{
			# replace missing values by NULL values for correct dispatch
			if( missing(method) ) method <- NULL
			if( missing(rank) ) rank <- NULL
			
			# apply KINOMO to the gene expression matrix			
			KINOMO(Biobase::exprs(x), rank, method, ...)
		}
	)
	
	

	setMethod('KINOMO', signature(x='matrix', rank='ExpressionSet', method='ANY'),
		function(x, rank, method, ...){
			# replace missing values by NULL values for correct dispatch
			if( missing(method) ) method <- NULL
			
			KINOMO(x, Biobase::exprs(rank), method, ...)
		}
	)
	
	
	

	setMethod('seed', signature(x='ExpressionSet', model='ANY', method='ANY'), 
		function(x, model, method, ...)
		{
			# replace missing values by NULL values for correct dispatch
			if( missing(method) ) method <- NULL
			if( missing(model) ) model <- NULL
			
			# apply KINOMO to the gene expression matrix			
			seed(Biobase::exprs(x), model, method, ...)
		}
	)
	
	#' Runs an KINOMO algorithm on the expression matrix of an \code{ExpressionSet} object.
	setMethod('run', signature(object='KINOMOStrategy', y='ExpressionSet', x='ANY'),
		function(object, y, x, ...){
			
			run(object, Biobase::exprs(y), x, ...)
			
		}
	)
		
	###% Method 'KINOMOModel' for 'ExpressionSet' target objects: 
	###% -> use the expression matrix of 'target' as the target matrix
	setMethod('KINOMOModel', signature(rank='ANY', target='ExpressionSet'),
			function(rank, target, ...){
				if( missing(rank) ) rank <- NULL
				# call KINOMOModel on the expression matrix
				KINOMOModel(rank, Biobase::exprs(target), ...)
			}	
	)
	setMethod('KINOMOModel', signature(rank='ExpressionSet', target='ANY'),
			function(rank, target, ...){
				if( missing(target) ) target <- NULL
				# call KINOMOModel on the expression matrix
				KINOMOModel(Biobase::exprs(rank), target, ...)
			}	
	)	
	
	###% Method 'rKINOMO' for 'ExpressionSet' target objects: 
	###% -> use the expression matrix of 'target' as the target matrix
	###% 
	setMethod('rKINOMO', signature(x='ANY', target='ExpressionSet'), 
		function(x, target, ...){
			rKINOMO(x, Biobase::exprs(target), ...)
		}
	)
	
	###% The method for an \code{ExpressionSet} object returns the data.frame that 
	###% contains the phenotypic data (i.e. \code{pData(object)})
	setMethod('.atrack', 'ExpressionSet', 
		function(object, data=NULL, ...){
			if( is.null(data) ) data <- t(Biobase::exprs(object))
			.atrack(Biobase::pData(object), data=data, ...)	
		}
	)
	
	

	setMethod('nneg', 'ExpressionSet'
			, function(object, ...){
				Biobase::exprs(object) <- nneg(Biobase::exprs(object), ...)
				object
			}
	)
	
	

	setMethod('rposneg', 'ExpressionSet'
			, function(object, ...){
				Biobase::exprs(object) <- rposneg(Biobase::exprs(object), ...)
				object
			}
	)
	
	###% Annotate the genes specific to each cluster.
	###%
	###% This function uses the \code{annaffy} package to generate an HTML table from the probe identifiers.


	## Assign BioConductor aliases
	###% number of metagenes
	nmeta <- nbasis
	###% get/set methods of basis matrix
	metagenes <- basis
	`metagenes<-` <- `basis<-`
	###% get/set methods of mixture coefficients matrix
	metaprofiles <- coef
	`metaprofiles<-` <- `coef<-`
	
	###% Get/Set methods for rows/columns names of the basis and mixture matrices
	# using the Biobase definition standard generics
	setGeneric('featureNames', package='Biobase')
	setGeneric('featureNames<-', package='Biobase')	
	setMethod('featureNames', 'KINOMO',
		function(object){
			rownames(object)
		}
	)
	setReplaceMethod('featureNames', 'KINOMO',
		function(object, value){
			rownames(object) <- value
			object
		}
	)
	###% For KINOMOfitX objects: returns the featureNames of the best fit 
	###% There is no replace method for KINOMOfitX objects
	setMethod('featureNames', 'KINOMOfitX',
		function(object){
			rownames(fit(object))
		}
	)
	
	setGeneric('sampleNames', package='Biobase')
	setGeneric('sampleNames<-', package='Biobase')	
	setMethod('sampleNames', 'KINOMO',
		function(object){
			colnames(object)
		}
	)
	setReplaceMethod('sampleNames', 'KINOMO',
		function(object, value){
			colnames(object) <- value
			object
		}
	)
	###% For KINOMOfitX objects: returns the sampleNames of the best fit 
	###% There is no replace method for KINOMOfitX objects
	setMethod('sampleNames', 'KINOMOfitX',
		function(object){
			colnames(fit(object))
		}
	)



	
	# return TRUE
	TRUE
}

}
