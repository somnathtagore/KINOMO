

#' @include options.R
#' @import grid
NULL


tryViewport <- function(name, verbose=FALSE){
	
	if( verbose ) message("vp - lookup for ", name)
	l <- grid.ls(viewports=TRUE, grobs=FALSE, print=FALSE)
	if( name %in% l$name ){
		downViewport(name)
	}	
}


current.vpPath_patched <- local({
    .current.vpPath <- NULL
    function(){
        
        f_current.vpPath <- .current.vpPath
        if( !.use.grid.patch() ) f_current.vpPath <- grid::current.vpPath
        else if( is.null(f_current.vpPath) ){ # load patch from installed file
            patch <- source(packagePath('scripts', 'grid.R', package = 'KINOMO'), local = TRUE)
            .current.vpPath <<- patch$value
            f_current.vpPath <- .current.vpPath
        }
        # call 
        f_current.vpPath()
    }
})

# Add new option to enable/disable grid patch
.OPTIONS$newOptions(grid.patch = FALSE)

#' \code{.use.grid.patch} tells if the user enabled patching grid.
#' @rdname grid
.use.grid.patch <- function(){
    !isCHECK() && KINOMO.getOption('grid.patch')   
}
