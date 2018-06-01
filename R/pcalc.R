

#' Calculation of the p-value for the test of H0: beta = 0 versus H1: beta > 0 for exact clustered logistic regression.  
#'
#' This function calculates the p-value for the test of H0: beta = 0 versus H1: beta > 0 for exact clustered logistic regression.  
#' @param zadmis The matrix of admissable vectors.  
#' @param nvec The vector of numbers of observations per cluster.  
#' @param z The vector of observed z-values for each cluster.  
#' @param x A vector of predictor values for each cluster.  
#' @keywords Exact Clustered Logistic Regression
#' @export
#' @examples
#' pcalc()


pcalc=function(zadmis,nvec,z,x){
	#calculate denominator
	denom=0
	nadmis = dim(zadmis)[2]

	for(i in 1:nadmis)
	{
		denom = denom + comprod(nvec, zadmis[,i])
	}

	#calculate numberator
	zextreme = extremez(zadmis,nvec,z,x)
	num=0
	nextreme = dim(zextreme)[2]

	for(i in 1:nextreme)
	{
		num = num + comprod(nvec,zextreme[,i])
	}

	p = num/denom

	p
} # end pcalc function
