
#' @import pkgmaker
#' @import registry
KINOMORegistry <- function(...) pkgmaker::packageRegistry(...)

# Returns the names of all the packages that contibute to all or a given
# package's primary registry  
registryContributors <- function(package, regname = NULL){
    regs <- packageRegistries(regname = regname, package = package, primary = TRUE)
    if( length(regs) ) unique(names(unlist(lapply(paste0(package, '::', regs), packageRegistries))))
}


KINOMOGet <- function(regname, name=NULL, ...){
	
	# retrieve from the given package's sub-registry
	pkgmaker::pkgreg_fetch(regname, key=name, ...)
	
}



###%
setGeneric('KINOMORegister', function(key, method, ...) standardGeneric('KINOMORegister') )
setMethod('KINOMORegister', signature(key='character'), 
	function(key, method, regname, ...){		
		#TODO: add functionality to save the registered strategy into a file for use is other R sessions
		
		parent.method <- attr(method, 'parent')
		tmpl <- if( !is.null(parent.method) && parent.method != key ){
			str_c(" based on template '", parent.method, "'")
		}
		setPackageRegistryEntry(regname, key, method, ..., where='KINOMO', msg=tmpl)
	}
)


