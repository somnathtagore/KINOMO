#' @include registry-algorithms.R
NULL




#' 
setGeneric('fcnnls', function(x, y, ...) standardGeneric('fcnnls') )

#' 
setMethod('fcnnls', signature(x='matrix', y='matrix'), 
	function(x, y, verbose=FALSE, pseudo=TRUE, ...){
		# load corpcor if necessary
		if( isTRUE(pseudo) ){ 
			library(corpcor)
		}
		
		# call the internal function
		res <- .fcnnls(x, y, verbose=verbose, pseudo=pseudo, ...)
		
		# process the result
		f <- x %*% res$coef
		resid <- y - f
        # set dimnames
        if( is.null(rownames(res$coef)) ) rownames(res$coef) <- colnames(x)

		# wrap up the result
		out <- list(x=res$coef, fitted=f, residuals=resid, deviance=norm(resid, 'F')^2, passive=res$Pset, pseudo=pseudo)
		class(out) <- 'fcnnls'
		out
	}
)
#' Shortcut for \code{fcnnls(as.matrix(x), y, ...)}.
setMethod('fcnnls', signature(x='numeric', y='matrix'), 
	function(x, y, ...){
		fcnnls(as.matrix(x), y, ...)
	}
)
#' Shortcut for \code{fcnnls(x, as.matrix(y), ...)}.
setMethod('fcnnls', signature(x='ANY', y='numeric'), 
	function(x, y, ...){
		fcnnls(x, as.matrix(y), ...)
	}
)

#' @export
print.fcnnls <- function(x, ...){
	cat("<object of class 'fcnnls': Fast Combinatorial Nonnegative Least Squares>\n")
	cat("Dimensions:", nrow(x$x)," x ", ncol(x$x), "\n")
	cat("Residual sum of squares:", x$deviance,"\n")
	cat("Active constraints:", length(x$passive)-sum(x$passive),"/", length(x$passive), "\n")
	cat("Inverse method:", 
			if( isTRUE(x$pseudo) ) 'pseudoinverse (corpcor)'
			else if( is.function(x$pseudo) ) str_fun(x$pseudo) 
			else 'QR (solve)', "\n")
	invisible(x)
}




.fcnnls <- function(x, y, verbose=FALSE, pseudo=FALSE, eps=0){
	
	# check arguments
	if( any(dim(y) == 0L) ){
		stop("Empty target matrix 'y' [", paste(dim(y), collapse=' x '), "]")
	}
	if( any(dim(x) == 0L) ){
		stop("Empty regression variable matrix 'x' [", paste(dim(x), collapse=' x '), "]")
	}
	
	# map arguments
	C <- x
	A <- y
# NNLS using normal equations and the fast combinatorial strategy
	#
	# I/O: [K, Pset] = fcnnls(C, A);
	# K = fcnnls(C, A);
	#
	# C is the nObs x lVar coefficient matrix
	# A is the nObs x pRHS matrix of observations
	# K is the lVar x pRHS solution matrix
	# Pset is the lVar x pRHS passive set logical array
	#
	# M. H. Van Benthem and M. R. Keenan
	# Sandia National Laboratories
	#
	# Pset: set of passive sets, one for each column
	# Fset: set of column indices for solutions that have not yet converged
	# Hset: set of column indices for currently infeasible solutions
	# Jset: working set of column indices for currently optimal solutions
	#
	# Check the input arguments for consistency and initializeerror(nargchk(2,2,nargin))
	nObs = nrow(C); lVar = ncol(C);
	if ( nrow(A)!= nObs ) stop('C and A have imcompatible sizes')
	pRHS = ncol(A);
	W = matrix(0, lVar, pRHS);
	iter=0; maxiter=3*lVar;
	# Precompute parts of pseudoinverse
	#CtC = t(C)%*%C; CtA = t(C)%*%A;
	CtC = crossprod(C); CtA = crossprod(C,A);
	
	# Obtain the initial feasible solution and corresponding passive set
	K = .cssls(CtC, CtA, pseudo=pseudo);
	Pset = K > 0;
	K[!Pset] = 0;
	D = K;
	# which columns of Pset do not have all entries TRUE?
	Fset = which( colSums(Pset) != lVar );
	#V+# Active set algorithm for NNLS main loop
	oitr=0; # HKim
	while ( length(Fset)>0 ) {
		
		oitr=oitr+1; if ( verbose && oitr > 5 ) cat(sprintf("%d ",oitr));# HKim
		
		#Vc# Solve for the passive variables (uses subroutine below)				
		K[,Fset] = .cssls(CtC, CtA[,Fset, drop=FALSE], Pset[,Fset, drop=FALSE], pseudo=pseudo);
		
		# Find any infeasible solutions
		# subset Fset on the columns that have at least one negative entry
		Hset = Fset[ colSums(K[,Fset, drop=FALSE] < eps) > 0 ];
		#V+# Make infeasible solutions feasible (standard NNLS inner loop)
		if ( length(Hset)>0 ){
			nHset = length(Hset);
			alpha = matrix(0, lVar, nHset);
			while ( nHset>0  && (iter < maxiter) ){
				iter = iter + 1; 
				alpha[,1:nHset] = Inf;
				#Vc# Find indices of negative variables in passive set
				ij = which( Pset[,Hset, drop=FALSE] & (K[,Hset, drop=FALSE] < eps) , arr.ind=TRUE);			
				i = ij[,1]; j = ij[,2]
				if ( length(i)==0 ) break;			
				hIdx = (j - 1) * lVar + i; # convert array indices to indexes relative to a lVar x nHset matrix
				negIdx = (Hset[j] - 1) * lVar + i; # convert array indices to index relative to the matrix K (i.e. same row index but col index is stored in Hset)
				
				alpha[hIdx] = D[negIdx] / (D[negIdx] - K[negIdx]);				
				alpha.inf <- alpha[,1:nHset, drop=FALSE]
				minIdx = max.col(-t(alpha.inf)) # get the indce of the min of each row
				alphaMin = alpha.inf[minIdx + (0:(nHset-1) * lVar)]
				alpha[,1:nHset] = matrix(alphaMin, lVar, nHset, byrow=TRUE);
				D[,Hset] = D[,Hset, drop=FALSE] - alpha[,1:nHset, drop=FALSE] * (D[,Hset, drop=FALSE]-K[,Hset, drop=FALSE]);			
				idx2zero = (Hset - 1) * lVar + minIdx; # convert array indices to index relative to the matrix D
				D[idx2zero] = 0;
				Pset[idx2zero] = FALSE;
				K[, Hset] = .cssls(CtC, CtA[,Hset, drop=FALSE], Pset[,Hset, drop=FALSE], pseudo=pseudo);
				# which column of K have at least one negative entry?
				Hset = which( colSums(K < eps) > 0 );
				nHset = length(Hset);
			}
		}
		#V-#
		
		#Vc# Make sure the solution has converged
		#if iter == maxiter, error('Maximum number iterations exceeded'), end
		# Check solutions for optimality
		W[,Fset] = CtA[,Fset, drop=FALSE] - CtC %*% K[,Fset, drop=FALSE];
		# which columns have all entries non-positive
		Jset = which( colSums( (ifelse(!(Pset[,Fset, drop=FALSE]),1,0) * W[,Fset, drop=FALSE]) > eps ) == 0 );
		Fset = setdiff(Fset, Fset[Jset]);
		
		if ( length(Fset) > 0 ){				
			#Vc# For non-optimal solutions, add the appropriate variable to Pset						
			# get indice of the maximum in each column
			mxidx = max.col( t(ifelse(!Pset[,Fset, drop=FALSE],1,0) * W[,Fset, drop=FALSE]) )
			Pset[ (Fset - 1) * lVar + mxidx ] = TRUE;
			D[,Fset] = K[,Fset, drop=FALSE];
		}		
	}
	#V-#
	
	# return K and Pset
	list(coef=K, Pset=Pset)
}
# ****************************** Subroutine****************************
#library(corpcor)
.cssls <- function(CtC, CtA, Pset=NULL, pseudo=FALSE){
	
	# use provided function
	if( is.function(pseudo) ){
		pseudoinverse <- pseudo
		pseudo <- TRUE
	}
	
	# Solve the set of equations CtA = CtC*K for the variables in set Pset
	# using the fast combinatorial approach
	K = matrix(0, nrow(CtA), ncol(CtA));	
	if ( is.null(Pset) || length(Pset)==0 || all(Pset) ){		
		K <- (if( !pseudo ) solve(CtC) else pseudoinverse(CtC)) %*% CtA;
		# K = pseudoinverse(CtC) %*% CtA;
		#K=pinv(CtC)*CtA;
	}else{
		lVar = nrow(Pset); pRHS = ncol(Pset);
		codedPset = as.numeric(2.^(seq(lVar-1,0,-1)) %*% Pset);
		sortedPset = sort(codedPset)
		sortedEset = order(codedPset)
		breaks = diff(sortedPset);
		breakIdx = c(0, which(breaks > 0 ), pRHS);
		for( k in seq(1,length(breakIdx)-1) ){
			cols2solve = sortedEset[ seq(breakIdx[k]+1, breakIdx[k+1])];
			vars = Pset[,sortedEset[breakIdx[k]+1]];			
			K[vars,cols2solve] <- (if( !pseudo ) solve(CtC[vars,vars, drop=FALSE]) else pseudoinverse(CtC[vars,vars, drop=FALSE])) %*% CtA[vars,cols2solve, drop=FALSE];
			#K[vars,cols2solve] <-  pseudoinverse(CtC[vars,vars, drop=FALSE])) %*% CtA[vars,cols2solve, drop=FALSE];
			#TODO: check if this is the right way or needs to be reversed
			#K(vars,cols2solve) = pinv(CtC(vars,vars))*CtA(vars,cols2solve);
		}
	}
	
	# return K
	K
}


#function [W,H,i] 
KINOMO_sKINOMO <- function(A, x, maxIter= KINOMO.getOption('maxIter') %||% 20000L, eta=-1, beta=0.01, bi_conv=c(0, 10), eps_conv=1e-4, version=c('R', 'L'), verbose=FALSE){
#KINOMOsh_comb <- function(A, k, param, verbose=FALSE, bi_conv=c(0, 10), eps_conv=1e-4, version=c('R', 'L')){
	
	# depending on the version: 
	# in version L: A is transposed while W and H are swapped and transposed
	version <- match.arg(version)
	if( version == 'L' ) A <- t(A) 
	#if( missing(param) ) param <- c(-1, 0.01)
	
	m = nrow(A); n = ncol(A); erravg1 = numeric();
	
	#eta=param[1]; beta=param[2]; 
	maxA=max(A); if ( eta<0 ) eta=maxA;
	eta2=eta^2;
	
	# bi_conv
	if( length(bi_conv) != 2 )
		stop("SKINOMO/", version, "::Invalid argument 'bi_conv' - value should be a 2-length numeric vector")
	wminchange=bi_conv[1]; iconv=bi_conv[2];
	
	## VALIDITY of parameters
	# eps_conv
	if( eps_conv <= 0 )
		stop("SKINOMO/", version, "::Invalid argument 'eps_conv' - value should be positive")
	# wminchange
	if( wminchange < 0 )
		stop("SKINOMO/", version, "::Invalid argument 'bi_conv' - bi_conv[1] (i.e 'wminchange') should be non-negative")
	# iconv
	if( iconv < 0 )
		stop("SKINOMO/", version, "::Invalid argument 'bi_conv' - bi_conv[2] (i.e 'iconv') should be non-negative")
	# beta
	if( beta <=0 )
		stop("SKINOMO/", version, "::Invalid argument 'beta' - value should be positive")
	##
	
	# initialize random W if no starting point is given
	if( isNumber(x) ){
		# rank is given by x
		k <- x
		message('# NOTE: Initialise W internally (runif)')
		W <- matrix(runif(m*k), m,k);	
		x <- NULL
	} else if( is.KINOMO(x) ){
		# rank is the number of basis components in x
		k <- nbasis(x)
		# seed the method (depends on the version to run)
		start <- if( version == 'R' ) basis(x) else t(coef(x))
		# check compatibility of the starting point with the target matrix
		if( any(dim(start) != c(m,k)) )
			stop("SKINOMO/", version, " - Invalid initialization - incompatible dimensions [expected: ", paste(c(m,k), collapse=' x '),", got: ", paste(dim(start), collapse=' x '), " ]")	
		# use the supplied starting point
		W <- start
	}else{
		stop("SKINOMO/", version, ' - Invalid argument `x`: must be a single numeric or an KINOMO model [', class(x), ']')
	}
	
	if ( verbose )
		cat(sprintf("--\nAlgorithm: SKINOMO/%s\nParameters: k=%d eta=%.4e beta (for sparse H)=%.4e wminchange=%d iconv=%d\n",
				version, k,eta,beta,wminchange,iconv));

	idxWold=rep(0, m); idxHold=rep(0, n); inc=0;
		
	# check validity of seed
	if( any(NAs <- is.na(W)) )
		stop("SKINOMO/", version, "::Invalid initialization - NAs found in the ", if(version=='R') 'basis (W)' else 'coefficient (H)' , " matrix [", sum(NAs), " NAs / ", length(NAs), " entries]")
	
	# normalize columns of W
	W= apply(W, 2, function(x) x / sqrt(sum(x^2)) );	

	I_k=diag(eta, k); betavec=rep(sqrt(beta), k); nrestart=0;
	i <- 0L
	while( i < maxIter){
		i <- i + 1L
		
		# min_h ||[[W; 1 ... 1]*H  - [A; 0 ... 0]||, s.t. H>=0, for given A and W.
	  	res = .fcnnls(rbind(W, betavec), rbind(A, rep(0, n)));	  	  	
		H = res[[1]]

		if ( any(rowSums(H)==0) ){
			if( verbose ) cat(sprintf("iter%d: 0 row in H eta=%.4e restart!\n",i,eta));
			nrestart=nrestart+1;
			if ( nrestart >= 10 ){
				warning("KINOMO::sKINOMO - Too many restarts due to too big 'beta' value [Computation stopped after the 9th restart]");
				break;
			}
						
			# re-initialize random W
			idxWold=rep(0, m); idxHold=rep(0, n); inc=0; 
			erravg1 <- numeric();# re-initialize base average error
			W=matrix(runif(m*k), m,k);
			W= apply(W, 2, function(x) x / sqrt(sum(x^2)) );  # normalize columns of W	
			next;
		}
		
		# min_w ||[H'; I_k]*W' - [A'; 0]||, s.t. W>=0, for given A and H. 
		res = .fcnnls(rbind(t(H), I_k), rbind(t(A), matrix(0, k,m))); 
		Wt = res[[1]]
		W= t(Wt);		

		# track the error (not computed unless tracking option is enabled in x)
		if( !is.null(x) ) 
			x <- trackError(x, .sKINOMO.objective(A, W, H, eta, beta), niter=i)
		
		# test convergence every 5 iterations OR if the base average error has not been computed yet
		if ( (i %% 5==0)  || (length(erravg1)==0) ){
			# indice of maximum for each row of W
			idxW = max.col(W)
			# indice of maximum for each column of H
			idxH = max.col(t(H))
			changedW=sum(idxW != idxWold); changedH=sum(idxH != idxHold);
			if ( (changedW<=wminchange) && (changedH==0) ) inc=inc+1
			else inc=0

			resmat=pmin(H, crossprod(W) %*% H - t(W) %*% A + matrix(beta, k , k) %*% H); resvec=as.numeric(resmat);
			resmat=pmin(W, W %*% tcrossprod(H) - A %*% t(H) + eta2 * W); resvec=c(resvec, as.numeric(resmat));
			conv=sum(abs(resvec)); #L1-norm      
			convnum=sum(abs(resvec)>0);
			erravg=conv/convnum;
			# compute base average error if necessary
			if ( length(erravg1)==0 )
				erravg1=erravg;
			
			if ( verbose && (i %% 1000==0) ){ # prints number of changing elements
				if( i==1000 ) cat("Track:\tIter\tInc\tchW\tchH\t---\terravg1\terravg\terravg/erravg1\n")
				cat(sprintf("\t%d\t%d\t%d\t%d\t---\terravg1: %.4e\terravg: %.4e\terravg/erravg1: %.4e\n",
					 i,inc,changedW,changedH,erravg1,erravg,erravg/erravg1));
			}
			
			#print(list(inc=inc, iconv=iconv, erravg=erravg, eps_conv=eps_conv, erravg1=erravg1))
			if ( (inc>=iconv) && (erravg<=eps_conv*erravg1) ) break;
			idxWold=idxW; idxHold=idxH; 
		}

	}
	
	if( verbose ) cat("--\n")
	
	# force to compute last error if not already done
	if( !is.null(x) ) 
		x <- trackError(x, .sKINOMO.objective(A, W, H, eta, beta), niter=i, force=TRUE)	

	# transpose and reswap the roles
	if( !is.null(x) ){ 
		if( version == 'L' ){
			.basis(x) <- t(H)
			.coef(x) <- t(W)
		}
		else{ 
			.basis(x) <- W 
			.coef(x) <- H
		}
		# set number of iterations performed
		niter(x) <- i
		
		return(x)	
	}else{
		res <- list(W=W, H=H)
		if( version == 'L' ){
			res$W <- t(H)
			res$H <- t(W)
		}		
		return(invisible(res))
	}
}


###% Computes the objective value for the SKINOMO algorithm
.sKINOMO.objective <- function(target, w, h, eta, beta){
	
	1/2 * ( sum( (target - (w %*% h))^2 ) 
				+ eta * sum(w^2) 
				+ beta * sum( colSums( h )^2 )
				)
}

sKINOMO.objective <- function(x, y, eta=-1, beta=0.01){
	.sKINOMO.objective(y, .basis(x), .coef(x), eta, beta)
}

###% Wrapper function to use the SKINOMO/R algorithm with the KINOMO package.
###%
.sKINOMO <- function(target, seed, maxIter=20000L, eta=-1, beta=0.01, bi_conv=c(0, 10), eps_conv=1e-4, ...){	
	
	# retrieve the version of SKINOMO algorithm from its name: 
	# it is defined by the last letter in the method's name (in upper case)
	name <- algorithm(seed)
	version <- toupper(substr(name, nchar(name), nchar(name)))
	
	# perform factorization using Kim and Park's algorithm
	ca <- match.call()
	ca[[1L]] <- as.name('KINOMO_sKINOMO')
	# target
	ca[['A']] <- ca[['target']]
	ca[['target']] <- NULL
	# seed
	ca[['x']] <- ca[['seed']]
	ca[['seed']] <- NULL
	# version 
	ca[['version']] <- version
	# verbose
	ca[['verbose']] <- verbose(seed)
	e <- parent.frame()
	sol <- eval(ca, envir=e)
#	KINOMO_sKINOMO(target, seed, ..., version = version, verbose = verbose(seed))
		
	# return solution
	return(sol)
}


KINOMOAlgorithm.SKINOMO_R <- setKINOMOMethod('sKINOMO/r', .sKINOMO, objective=sKINOMO.objective)
#' @aliases SKINOMO/L-KINOMO
#' @rdname SKINOMO-KINOMO
KINOMOAlgorithm.SKINOMO_L <- setKINOMOMethod('sKINOMO/l', .sKINOMO, objective=sKINOMO.objective)
