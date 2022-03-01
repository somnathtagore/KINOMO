#' @import graphics
#' @import rngtools
#' @import digest
#' @import stringr
#' @import stats
#' @import methods
NULL
#library(digest)

#' Defunct Functions and Classes in the KINOMO Package
#' 
#' @name KINOMO-defunct
#' @rdname KINOMO-defunct
NULL

#' Deprecated Functions in the Package KINOMO
#' 
#' @param object an R object
#' @param ... extra arguments 
#' 
#' @name KINOMO-deprecated
#' @rdname KINOMO-deprecated
NULL




#' 
NA

devKINOMO <- function(){
	.LOCAL_PKG_NAME <- 'KINOMO'
	requireNamespace('devtools')
	devtools::load_all(.LOCAL_PKG_NAME)
	compile_src(.LOCAL_PKG_NAME)
}

# local config info
KINOMOConfig <- mkoptions()

.onLoad <- function(libname, pkgname) {
		
	# set default number of cores
	if( pkgmaker::isCHECK() ){
		options(cores=2)
	}else{
		if( nchar(nc <- Sys.getenv('R_PACKAGE_KINOMO_CORES')) > 0 ){
			try({
				KINOMO.options(cores=as.numeric(nc))
			})
		}   
	}
    # use grid patch?
    KINOMO.options(grid.patch = !isFALSE(Sys.getenv_value('R_PACKAGE_KINOMO_GRID_PATCH')))
    
    pkgEnv <- pkgmaker::packageEnv()
	.init.sequence <- function(){
	
		## 0. INITIALIZE PACKAGE SPECFIC OPTIONS
		#.init.KINOMO.options()
				
		## 1. INITIALIZE THE KINOMO MODELS
		.init.KINOMO.models()		
		
		## 2. INITIALIZE BIOC LAYER
		b <- body(.onLoad.KINOMO.bioc)
		bioc.loaded <- eval(b, envir=pkgEnv)
		KINOMOConfig(bioc=bioc.loaded)
		
		# 3. SHARED MEMORY
		if( .Platform$OS.type != 'windows' ){
			msg <- if( !require.quiet('bigmemory', character.only=TRUE) ) 'bigmemory'
					else if( !require.quiet('synchronicity', character.only=TRUE) ) 'synchronicity'
					else TRUE
			
			KINOMOConfig(shared.memory=msg)
		}
		#
	}
		
	# run intialization sequence suppressing messages or not depending on verbosity options
	.init.sequence()
	if( getOption('verbose') ) .init.sequence()
	else suppressMessages(.init.sequence())
	
	
	return(invisible())
}

.onUnload <- function(libpath) {
	
	# unload compiled library
	dlls <- names(base::getLoadedDLLs())
	if ( 'KINOMO' %in%  dlls )
		library.dynam.unload("KINOMO", libpath);	
}

.onAttach <- function(libname, pkgname){
	
	# build startup message
	msg <- NULL
	details <- NULL
	## 1. CHECK BIOC LAYER
	bioc.loaded <- KINOMOConfig('bioc')[[1L]]
	msg <- paste0(msg, 'BioConductor layer')
	if( is(bioc.loaded, 'try-error') ) msg <- paste0(msg, ' [ERROR]')
	else if ( isTRUE(bioc.loaded) ) msg <- paste0(msg, ' [OK]')
	else{
		msg <- paste0(msg, ' [NO: missing Biobase]')
		details <- c(details, "  To enable the Bioconductor layer, try: install.extras('", pkgname, "') [with Bioconductor repository enabled]")
	}
	
	# 2. SHARED MEMORY
	msg <- paste0(msg, ' | Shared memory capabilities')
	if( .Platform$OS.type != 'windows' ){
		conf <- KINOMOConfig('shared.memory')[[1L]]
		if( isTRUE(conf) ) msg <- paste0(msg, ' [OK]')
		else{
			msg <- paste0(msg, ' [NO: ', conf, ']')
			details <- c(details, "  To enable shared memory capabilities, try: install.extras('", pkgname, "')")
		}
	}else msg <- paste0(msg, ' [NO: windows]')
	#
	
	# 3. NUMBER OF CORES
	msg <- paste0(msg, ' | Cores ', getMaxCores(), '/', getMaxCores(limit=FALSE))
	#
	
	# FINAL. CRAN FLAG
	if( pkgmaker::isCHECK() ){
		msg <- paste0(msg, ' | CRAN check')
	}
	#
	
	# print startup message
	ver <- if( isDevNamespace() ){
		paste0(' [', utils::packageVersion(pkgname), '-devel', ']') 
	}#else{
#		utils::packageVersion(pkgname, lib.loc = libname)
#	}
	packageStartupMessage(pkgname, ver, ' - ', msg)
	if( !is.null(details) ){
		packageStartupMessage(paste(details, collapse="\n"))
	}
}

