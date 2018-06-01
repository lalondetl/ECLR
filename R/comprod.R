

#' Calculation of the product of combinations used in conditional probability distribution for exact clustered logistic regression.  
#'
#' This function calculates the product of combinations used in conditional probability distribution for exact clustered logistic regression.  
#' @param nvec The vector of numbers of observations per cluster.  
#' @param zvec The vector of observed z-values (the sum of successses in cluster)
#' @keywords Exact Clustered Logistic Regression
#' @export
#' @examples
#' comprod()


comprod = function(nvec,zvec){
	prod=1
	N = length(nvec)

	for(i in 1:N)
	{
		prod = prod*(choose(nvec[i],zvec[i]))
	}

	prod
} # end comprod function
