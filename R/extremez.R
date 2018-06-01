

#' Calculation of all admissable z-vectors with corresponding t-statistics at least as large as the observed value of t.  
#'
#' This function performs the calculation of all admissable z-vectors with corresponding t-statistics at least as large as the observed value of t for exact clustered logistic regression.  
#' @param zadmis The matrix of admissable vectors.  
#' @param nvec The vector of numbers of observations per cluster.  
#' @param z The vector of observed z-values.  
#' @param x A vector of predictor values for each cluster.  
#' @keywords Exact Clustered Logistic Regression
#' @export
#' @examples
#' extremez()


extremez=function(zadmis,nvec,z,x){
	ncols = dim(zadmis)[2]
	tobs = tcalc(z,x)

	#initialize answer to include observed values
	zextreme = z

	for(i in 1:ncols) #run through all admissable vectors
	{
		zcurrent = zadmis[,i]
		tcurrent = tcalc(zcurrent,x)

		if(tcurrent >= tobs)
		{
			zextreme = cbind(zextreme,zcurrent)
		}
	} #end for

	#remove first column; observed z-values are recorded twice
	zextremefinal = zextreme[,2:dim(zextreme)[2]]

	zextremefinal
} # end extremez function
