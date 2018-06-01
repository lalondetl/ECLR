

#' Initialization of admissibility test recursion.  
#'
#' This function initializes the admissability test recursion for exact clustered logistic regression.  The matrix of admissable vectors starts as observed values, zvec.  The final output is a matrix whose columns are z-vectors that give the same sufficient statstics as those observed.  
#' @param nvec The vector of numbers of observations per cluster.  
#' @param zvec The vector of observed z-values (the sum of successses in cluster)
#' @keywords Exact Clustered Logistic Regression
#' @export
#' @examples
#' admis()


admis = function(nvec,zvec){
	#N is the number of clusters
	N = length(nvec)
	
	#s1 and s2 are the observed values of the sufficient statistics
	s1 = s1calc(zvec)
	s2 = s2calc(zvec,nvec)

	#begin with first cluster
	k = 1
	n1 = nvec[1]

	#initially zmat is the observed vector of values
	zmat = zvec

	#run through all possible numbers of successes in first cluster
	for(z in 0:n1)
	{
		znew=as.vector(z)
		
		if(N==1) #check case of one cluster
		{
			s1current = s1calc(znew)
			s2current = s2calc(znew,nvec)

			if(s1current==s1 && s2current==s2)
			{
				zmat = cbind(zmat,znew)
			}

		} #end one cluster check

		else #more than one cluster
		{
			#apply recursive function to find full admissable vectors
			zmat = admisrec(nvec,zmat,znew,s1,s2)
		} #zmat will have all admissable vectors beginning with znew added as columns

	} #end for

	#zmat will have all admissable vectors for all values of z added as columns
	#observed z-vector will occur twice due to starting the process with this vector
	#zmatfinal has this first vector removed
	ncols=dim(zmat)[2]
	zmatfinal=zmat[,2:ncols]

	zmatfinal

} # end admis function