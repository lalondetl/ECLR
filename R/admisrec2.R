

#' Admissibility test recursion, two-level model.  
#'
#' This function performs the admissability test recursion for exact clustered logistic regression, two-level model.  The matrix of admissable vectors starts as observed values, zvec.  The final output is a matrix whose columns are z-vectors that give the same sufficient statstics as those observed.  
#' @param nvec The vector of numbers of observations per cluster.  
#' @param zmat The vector of observed z-values.  
#' @param zvec The vector of (possibly incomplete) z-values.  
#' @param jvec The numbers of inner clusters in each outer cluster.  
#' @param s1 Observed value of the first sufficient statistic.  
#' @param s2 Observed value of the second sufficient statistic.  
#' @param s2 Observed value of the third sufficient statistic.  
#' @keywords Exact Clustered Logistic Regression
#' @export
#' @examples
#' admisrec2()


admisrec2 = function(nvec,zmat,zvec,jvec,s1,s2,s3){
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
			zmatfinal = admisrec2(nvec,zmatfinal,zvecnew,jvec,s1,s2,s3)
		}

		else #length=N, zvecnew is a full vector of potential observations
		{
			s1current = s1calc(zvecnew)
			s2current = s2calc(zvecnew,nvec)
			s3current = s3calc(zvecnew,nvec,jvec)

			if(s1current==s1 && s2current==s2 && s3current==s3)
			{
				zmatfinal = cbind(zmatfinal,zvecnew)
			}
		} #end else
		#zmatfinal contains all admissable vectors with z in component i
	} #end for
	
	#zmatfinal contains all admissable vectors for all values of z in component i
	zmatfinal

} # end admisrec2 function
