

#' Calculation of third sufficient statistic for exact clustered logistic regression.  
#'
#' This function calculates the third sufficient statistic, s3, for exact clustered logistic regression.
#' @param zvec The vector of observed z-values (the sum of successses in cluster)
#' @param nvec The vector of numbers of observations per cluster
#' @param jvec The vector of numbers of inner clusters per outer cluster
#' @keywords Exact Clustered Logistic Regression
#' @export
#' @examples
#' s3calc()


s3calc = function(zvec,nvec,jvec){
	k = length(jvec) #number of outer clusters
	s3 = 0

	count = 1

	for(i in 1:k) #run through all inner
	{
		nsum = 0
		zsum = 0

		for(j in count:(count+jvec[i]-1))
		{
			nsum = nsum + nvec[j]
			zsum = zsum + zvec[j]
		}

		s3i = 0

		for(j in count:(count+jvec[i]-1))
		{
			s3i = s3i + zvec[j]*(nsum - zsum)
		}

		s3 = s3 + s3i
		count = count + jvec[i]
	} #end for

	s3
} # end s3calc function