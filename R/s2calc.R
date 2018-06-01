

#' Calculation of second sufficient statistic for exact clustered logistic regression.  
#'
#' This function calculates the second sufficient statistic, s2, for exact clustered logistic regression.
#' @param zvec The vector of observed z-values (the sum of successses in cluster)
#' @param nvec The vector of numbers of observations per cluster
#' @keywords Exact Clustered Logistic Regression
#' @export
#' @examples
#' s2calc()


s2calc = function(zvec, nvec){
	s2 = sum(zvec*(nvec-zvec))

	s2
} #end s2calc