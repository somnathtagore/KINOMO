
#' @include options.R
#' @import grid
NULL

#' Internal Grid Extension
#' 
#' These functions enable mixing base and grid graphics in \code{\link{aheatmap}},
#' by avoiding calls to the grid internal function \code{'L_gridDirty'}.
#' They are not exported (i.e. not tampering core functions) and are only meant for internal 
#' use within the \pkg{NMF} package.
#'
#' \code{tryViewport} tries to go down to a viewport in the current tree, 
#' given its name.
#' 
#' @details
#' \code{tryViewport} uses \code{\link[grid]{grid.ls}} and not 
#' \code{\link{seekViewport}} as the latter would reset the graphic device 
#' and break the mix grid/base graphic capability.  
#' 
#' @param name viewport name 
#' @param verbose toggle verbosity
#' 
#' @rdname grid
#' @keywords internal
tryViewport <- function(name, verbose=FALSE){
	
	if( verbose ) message("vp - lookup for ", name)
	l <- grid.ls(viewports=TRUE, grobs=FALSE, print=FALSE)
	if( name %in% l$name ){
		downViewport(name)
	}	
}

#' \code{current.vpPath_patched} aims at substituting \code{\link[grid]{current.vpPath}}, 
#' so that the graphic engine is not reset. 
#' This is essentially to prevent outputting a blank page at the beginning of PDF 
#' graphic engines. 
#' 
#' @rdname grid
current.vpPath_patched <- local({
    .current.vpPath <- NULL
    function(force.unpatched = FALSE, verbose = FALSE){
        
        f_current.vpPath <- .current.vpPath
        if( force.unpatched || !.use.grid.patch() ) f_current.vpPath <- grid::current.vpPath
        else if( is.null(f_current.vpPath) ){ # load patch from installed file
            patch <- source(packagePath('scripts', 'grid.R', package = 'NMF'), local = TRUE)
            .current.vpPath <<- patch$value
            f_current.vpPath <- .current.vpPath
        }
        # call 
        if( verbose ){
          version <- c('patched', 'grid')[(digest(f_current.vpPath) == digest(grid::current.vpPath)) + 1L]
          message("Using current.vpPath version: ", version)
        }
        f_current.vpPath()
    }
})

#' @section Plotting:
#' \describe{
#' \item{grid.patch}{logical that indicates if annotated heatmaps should use a set of functions 
#' from the `grid` packages that were patched to enable mixing (`TRUE`) or stock `grid` functions 
#' (`FALSE`).
#' The default value for this option is `NULL`, which avoids the saving of heatmaps in pdf files 
#' to generate a blank first page}
#' }
#' 
#' @rdname options
#' @name options
NULL
# Add new option to enable/disable grid patch
.OPTIONS$newOptions(grid.patch = NULL)

#' \code{.use.grid.patch} tells if the user enabled patching grid.
#' @rdname grid
.use.grid.patch <- function(){
  # default is to not use patched grid unless within a pdf device
  default <- names(dev.cur()) %in% 'pdf'
  !isCHECK() && nmf.getOption('grid.patch', default)
}
