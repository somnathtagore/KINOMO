

#' @include KINOMO-class.R
NULL

if( FALSE ){ #START_DEACTIVATE
	
## #' Factory Method for KINOMO Model Classes
## #' 
## #' Defines two S4 classes for representing KINOMO models: one to hold data from 
## #' the actual model, the other one to hold fitting data for model estimated with 
## #' the function \code{\link{KINOMO}}. 
## #' 
#setKINOMOClass <- function(Class, ..., where=topns(), contains='KINOMOstd', VERBOSE=TRUE){
#	
#	# add 'KINOMO' to contains if necessary
#	wKINOMO <- sapply(contains, isKINOMOclass)
#	if( !length(wKINOMO) || !any(wKINOMO) ){
#		contains <- c(contains, 'KINOMOstd')
#		parentKINOMOClass <- 'KINOMOstd'
#	}else{
#		parentKINOMOClass <- contains[which(wKINOMO)]
#	}
#	
#	# extract KINOMO prefix if present
#	Class <- sub('^KINOMO(.*)', "\\1", Class)
#	# define class names
#	KINOMOClass <- str_c('KINOMO', Class)
#	KINOMOfitClass <- str_c(KINOMOClass, '_fit')
#	if( VERBOSE ){
#		message("Defining KINOMO classes: ", KINOMOClass , "(", parentKINOMOClass , ") and "
#				, KINOMOfitClass, ' in ', str_ns(where), ' ... '
#				, appendLF=FALSE)
#	}
#	# 1. Create model class
#	setClass(KINOMOClass, ..., where=where, contains=contains)
#	
#	# 2. Create model fit class (in the same environment as the model class)
#	e <- packageEnv(getClassDef(KINOMOClass)@package)
#	setClass(KINOMOfitClass, where=e, contains=c('KINOMOfit', KINOMOClass))
#	
#	if( VERBOSE ) message('OK')
#	
#	# return the name of the two new classes
#	c(KINOMOClass, KINOMOfitClass)
#}  

}#END_DEACTIVATE
