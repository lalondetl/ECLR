

#' Calculation of a sufficient statistic for a parameter of iterest for exact clustered logistic regression.  
#'
#' This function calculates a sufficient statistic, t, for a parameter of iterest for exact clustered logistic regression.
#' @param zvec The vector of observed z-values (the sum of successses in cluster)
#' @param x A vector of predictor values for each cluster.  
#' @keywords Exact Clustered Logistic Regression
#' @export
#' @examples
#' tcalc()


tcalc = function(zvec,x){
	t = sum(zvec*x)

	t
}