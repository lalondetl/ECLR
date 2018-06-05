/*****************************************************************/
/**	SAS Code: %exactlogistic 								    **/
/** Programmer: Kyle Irimata									**/
/** Description: Performs exact logistic regresion using		**/
/**		method in Troxler, Lalonde, and Wilson (2011)  			**/
/*****************************************************************/

%macro exactlogistic(data=, levels=1 );

proc iml;
*nvec - vector of observations per cluster;
*zmat - vector of observed z values;
*zvec - a possibly incomplete vector of z-values;
*s1, s2 - observed valuse of sufficient statistics;

*Read data into iml;
use &data;
read all var _all_ into dat;
close &data;

*Calculate the cluster sizes;
call tabulate(group,jcount,dat[,1]);

*s1calc;
start s1calc(zvec);
	s1 = sum(zvec);
	return (s1);
finish s1calc;

*s2calc;
start s2calc(zvec,nvec);
	s2 = sum(zvec`*(nvec-zvec));
	return (s2);
finish s2calc;

*s3calc;
*zvec is vector of observed z values (sum of successes in inner cluster);
*nvec, number of observations per inner cluster;
*jvec is number of inner clusters per outer cluster;
start s3calc(zvec,nvec,jvec);
	*number of outer clusters);
	k = ncol(jvec);	
	s3 = 0;

	count = 1;

	*Loop through all inner;
	do i=1 to k;
		nsum = 0;
		zsum = 0;
		
		endcount = count + jvec[i] - 1;
		do j=count to endcount;
			nsum = nsum + nvec[j];
			zsum = zsum + zvec[j];
		end;

		s3i = 0;

		do j=count to endcount;
			s3i = s3i + zvec[j]*(nsum - zsum);
		end;

		s3 = s3 + s3i;
		count = count + jvec[i];
	end;
	return (s3);
finish s3calc;

*tcalc;
start tcalc(zvec,x);
	t = sum(zvec`*x);
	return (t);
finish tcalc;

*Initializes the admissibility recursion for one level of clustering;
start admis(nvec, zvec);
	*Number of inner clusters;
	n = nrow(nvec);

	*obtain the observed sufficient statistics;
	s1 = s1calc(zvec);
	s2 = s2calc(zvec,nvec);

	*Begin with first cluster;
	k = 1;
	n1 = nvec[1];
	
	*Initially zmat is the observed vector of values;
	zmat = zvec;

	*Go through all possible numbers of successes in the first cluster;
	do z=0 to n1;
		znew = z;

		*Check for one cluster;
		if n=1 then do;
			s1current = s1calc(znew);
			s2current = s2calc(znew,nvec);

			if s1current = s1 & s2current = s2 then
				zmat= zmat||znew;
		end;
		*end check for one cluster;
		
		*for more than one cluster;
		else do;
			*Apply the recursion;
			*If there are errors, may need to define the admisrec, then assign it to zmat separately;
			zmat = admisrec(nvec, zmat, znew, s1, s2);
		end;
	end;
	
	cols = ncol(zmat);
	zmatfinal = zmat[,2:cols];

	return(zmatfinal);

finish admis;


*Recursive function to go through all admissable vectors;
*nvec is vector of observationsper cluster;
*zmat is vector of observed z-values;
*zvec a vector of z-values;
*s1,s2 are observed values of sufficient statistics;
start admisrec(nvec, zmat, zvec, s1, s2);
	n = nrow(nvec);

	zmatfinal = zmat;

	*Set the possible newvector to the initialized value;
	zvecnew = zvec;

	*continue until there is a vector of complete observations;
	do i=2 to n;
		temp = 0:nvec[i];

		templength = ncol(temp);
		zlength = ncol(zvecnew);


		zvecnew = repeat(zvecnew,1,templength);

		zreplength = ncol(zvecnew);


		temprep = j(1,zreplength,.);

		startind=1;
		do j=0 to nvec[i];
			endind=startind+zlength-1;
			temprep[,startind:endind] = j;
			startind = endind + 1;
		end;
		zvecnew = zvecnew // temprep;
	
	end;



	*Add the admissable columns to zmatfinal;
	allpos = ncol(zvecnew);
	do z=1 to allpos;
		s1current = s1calc(zvecnew[,z]);
		s2current = s2calc(zvecnew[,z],nvec);
		if s1current = s1 & s2current = s2 then
			zmatfinal= zmatfinal||zvecnew[,z];
	end;

	return(zmatfinal);
finish admisrec;



*Initializes the admissibility  recursion for two levels of clustering;
start admis2(nvec, zvec,jvec);
	*Number of inner clusters;
	n = nrow(nvec);

	*obtain the observed sufficient statistics;
	s1 = s1calc(zvec);
	s2 = s2calc(zvec,nvec);
	s3 = s3calc(zvec,nvec,jvec);

	*Begin with first cluster;
	k = 1;
	n1 = nvec[1];
	
	*Initially zmat is the observed vector of values;
	zmat = zvec;

	*Go through all possible numbers of successes in the first cluster;
	do z=0 to n1;
		znew = z;

		*Check for one cluster;
		if n=1 then do;
			s1current = s1calc(znew);
			s2current = s2calc(znew,nvec);
			s3current = s3calc(znew,nvec,jvec);

			if s1current = s1 & s2current = s2 & s3current = s3 then
				zmat= zmat||znew;
		end;
		*end check for one cluster;
		
		*for more than one cluster;
		else do;
			*Apply the recursion;
			*If there are errors, may need to define the admisrec, then assign it to zmat separately;
			zmat = admisrec2(nvec, zmat, znew,jvec, s1, s2, s3);
		end;
	end;
	
	cols = ncol(zmat);
	zmatfinal = zmat[,2:cols];

	return(zmatfinal);

finish admis2;


*Recursive function to go through all admissable vectors;
*nvec is vector of observationsper cluster;
*zmat is vector of observed z-values;
*zvec a vector of z-values;
*s1,s2 are observed values of sufficient statistics;
start admisrec2(nvec, zmat, zvec, jvec, s1, s2, s3);
	n = nrow(nvec);

	zmatfinal = zmat;

	*Set the possible newvector to the initialized value;
	zvecnew = zvec;

	*continue until there is a vector of complete observations;
	do i=2 to n;
		temp = 0:nvec[i];

		templength = ncol(temp);
		zlength = ncol(zvecnew);


		zvecnew = repeat(zvecnew,1,templength);

		zreplength = ncol(zvecnew);


		temprep = j(1,zreplength,.);

		startind=1;
		do j=0 to nvec[i];
			endind=startind+zlength-1;
			temprep[,startind:endind] = j;
			startind = endind + 1;
		end;
		zvecnew = zvecnew // temprep;

		*Check if the columns are possible under s1;
		s1check = zvecnew[+,];
		feascheck = loc(s1check <= s1);

		zvecnew = zvecnew[,feascheck];
	end;



	*Add the admissable columns to zmatfinal;
	allpos = ncol(zvecnew);
	do z=1 to allpos;
		s1current = s1calc(zvecnew[,z]);
		s2current = s2calc(zvecnew[,z],nvec);
		s3current = s3calc(zvecnew[,z],nvec,jvec);
		if s1current = s1 & s2current = s2 & s3current=s3 then
			zmatfinal= zmatfinal||zvecnew[,z];
	end;

	return(zmatfinal);
finish admisrec2;

*Calculates all admissible z-vectors with corresponding t-statistics at least as large as the observed t;
start extremez(zadmis,nvec,z,x);
	cols = ncol(zadmis);
	tobs = tcalc(z,x);

	*Initialize output to include observed values;
	zextreme = z;

	*Go through all admissible vectors;
	do i=1 to cols;
		zcurrent = zadmis[,i];
		tcurrent = tcalc(zcurrent,x);

		if(tcurrent >= tobs) then
			zextreme = zextreme||zcurrent;
	end;

	*Remove first column of observed z-values;
	nz = ncol(zextreme);
	zextremefinal = zextreme[,2:nz];

	return(zextremefinal);

finish extremez;

*Calculates the products of combinations used in the conditional probability distribution;
start comprod(nvec,zvec);
	prod=1;
	n = nrow(nvec);

	do i=1 to n;
		prod = prod*comb(nvec[i],zvec[i]);
	end;
	return(prod);
finish comprod;

*Calculates the p-value for the hypothesis test;
*H0: beta=0;
*Ha: beta>0;
start pcalc(zadmis,nvec,z,x);
	denom=0;
	nadmis = ncol(zadmis);

	do i=1 to nadmis;
		denom = denom + comprod(nvec,zadmis[,i]);
	end;
	
	zextreme = extremez(zadmis,nvec,z,x);
	num=0;
	nextreme = ncol(zextreme);

	do i=1 to nextreme;
		num = num + comprod(nvec,zextreme[,i]);
	end;

	p = num/denom;

	return(p);

finish pcalc;

if &levels=1 then do;
	admistest = admis(dat[,2],dat[,3]);
	pval = pcalc(admistest,dat[,2],dat[,3],dat[,4]);
end;

else do;
	admistest = admis2(dat[,2],dat[,3],jcount);
	pval = pcalc(admistest,dat[,2],dat[,3],dat[,4]);
end;


print pval;


quit;

%mend exactlogistic;
