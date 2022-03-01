

KINOMOReport <- function(x, rank, method, colClass = NULL, ..., output = NULL, template = NULL){
	
	requireNamespace('knitr')
	#library(knitr)
	if( is.null(template) )
		template <- system.file('scripts/report.Rmd', package = 'KINOMO')
	x <- force(x)
	rank <- force(rank)
	method <- force(method)
    if( isString(method) ) method <- list(method)
    
	args <- list(...)
	KINOMORun <- function(x, rank, method, ...){
		args <- expand_dots(args)
		str(args)
		do.call(KINOMO, c(list(x, rank, method), args))
	}
    accuracy <- NA
    res <- NA
	knitr::knit2html(template)
	res <- list(fits = res, accuracy = accuracy)
	saveRDS(res, file = 'report_results.rds')
	invisible(res)
}
