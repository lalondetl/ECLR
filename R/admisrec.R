

#' Admissibility test recursion.  
#'
#' This function performs the admissability test recursion for exact clustered logistic regression.  The matrix of admissable vectors starts as observed values, zvec.  The final output is a matrix whose columns are z-vectors that give the same sufficient statstics as those observed.  
#' @param nvec The vector of numbers of observations per cluster.  
#' @param zmat The vector of observed z-values.  
#' @param zvec The vector of (possibly incomplete) z-values.  
#' @param s1 Observed value of the first sufficient statistic.  
#' @param s2 Observed value of the second sufficient statistic.  
#' @keywords Exact Clustered Logistic Regression
#' @export
#' @examples
#' admisrec()



admisrec = function(nvec,zmat,zvec,s1,s2){
	N = length(nvec)

	#k is the current cluster to be analyzed
	k = length(zvec) + 1
	ni = nvec[k]

	zmatfinal = zmat #will be updated throughout function

	#run through all possible numbers of successes in cluster i
	for(z in 0:ni)
	{
		#add new z-value to vector
		zvecnew = rbind(zvec,z)

		#if this is not a full vector of potential observations, call again
		if(length(zvecnew) < N)
		{
			zmatfinal = admisrec(nvec,zmatfinal,zvecnew,s1,s2)
		}

		else #length=N, zvecnew is a full vector of potential observations
		{
			s1current = s1calc(zvecnew)
			s2current = s2calc(zvecnew,nvec)

			if(s1current==s1 && s2current==s2)
			{
				zmatfinal = cbind(zmatfinal,zvecnew)
			}
		} #end else
		#zmatfinal contains all admissable vectors with z in component i
	} #end for
	
	#zmatfinal contains all admissable vectors for all values of z in component i
	zmatfinal

} # end admisrec function
