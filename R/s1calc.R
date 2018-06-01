

#' Calculation of first sufficient statistic for exact clustered logistic regression.  
#'
#' This function calculates the first sufficient statistic, s1, for exact clustered logistic regression.
#' @param zvec The vector of observed z-values (the sum of successses in cluster)
#' @keywords Exact Clustered Logistic Regression
#' @export
#' @examples
#' s1calc()


s1calc = function(zvec){
	s1 = sum(zvec)

	s1
} #end s1calc