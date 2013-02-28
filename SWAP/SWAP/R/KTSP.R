SWAP.rank<-function(x)
{
	return((apply(x,2,rank)));

}



SWAP.CalculateSignedScore<-function(  situation , data1, data2 )
{
	
	n = length(situation);
	m1 = nrow(data1);
	m2 = nrow(data2);
	
	datatied = SWAP.rank(rbind(data1,data2));
	data1tied = datatied[1:m1,];
	data2tied = datatied[(m1+1):(m1+m2),];


	situationV = as.vector(situation);
	if( length(unique(situationV)) != 2)
		stop("The situation must contain only exact two different variable.")

	labels= unique(situationV);
	
	
	
	d<-
	.C(
	"CalculateSignedScoreCore",
	as.integer(as.numeric(situationV==labels[1])), as.integer(n),
	as.double(data1tied),as.integer(m1),
	as.double(data2tied),as.integer(m2),
	as.double(matrix(0,m1,m2)),
	as.double(matrix(0,m1,m2)),
	as.double(1),
	as.double(matrix(0,m1,m2)),
	as.double(1),
	as.double(matrix(0,m1,m2)),
	as.double(matrix(0,m1,m2))
	);
	
	score=(matrix(d[[7]],nrow=m1))/2;
#	a=(matrix(d[[7]],nrow=m1));
#     M=(matrix(d[[8]],nrow=1));
#	k=(matrix(d[[9]],nrow=m1));
#	N=(matrix(d[[10]],nrow=1));
     	P=(matrix(d[[12]],nrow=m1));
	Q=(matrix(d[[13]],nrow=m1));
	
	if( length(rownames(data1))>0 &&  length(rownames(data2))>0 )
	{
		names1 = rownames(data1);
		names2 = rownames(data2);
		
		rownames(score)<-names1;
		colnames(score)<-names2;
		
		rownames(P)<-names1;
		colnames(P)<-names2;
		
		rownames(Q)<-names1;
		colnames(Q)<-names2;	

	}
	retVal<-list(score=score,P=P,Q=Q,labels=labels);	

}


SWAP.CalculateSignedScore.Restricted<-function(  situation , data1, data2, pairs )
{
	
	n = length(situation);
	m1 = nrow(data1);
	m2 = nrow(data2);
	pairsno = nrow(pairs);	

	datatied = SWAP.rank(rbind(data1,data2));
	data1tied = datatied[1:m1,];
	data2tied = datatied[(m1+1):(m1+m2),];


	situationV = as.vector(situation);
	if( length(unique(situationV)) != 2)
		stop("The situation must contain only exact two different variable.")

	labels= unique(situationV);
	

	
	if( class(pairs[1,1])=="character" )
	{
		pairsind = matrix(0,pairsno ,2);
		pairsind[,1] = match(pairs[,1],rownames(data1));
		pairsind[,2] = match(pairs[,2],rownames(data2));
		if( sum(is.na( 	pairsind ))>0)
			stop("The genes in pairs do not show in data1 or data2.");	
	}
	else
	{
		
		if( min(pairs)<1 || max(pairs[,1])>nrow(data1) || max(pairs[,2])>nrow(data2))
			stop("The genes in pairs do not show in data1 or data2.");	
		pairsind = pairs;
	}
	
	d<-
	.C(
	"CalculateSignedScoreRestrictedPairsCore",
	as.integer(as.numeric(situationV==labels[1])),
	as.integer(n),
	as.double(data1tied),as.integer(m1),
	as.double(data2tied),as.integer(m2),
	as.integer(pairsind[,1]-1),
	as.integer(pairsind [,2]-1),
	as.integer(pairsno),
	as.double(matrix(0,pairsno,1)),
	as.double(matrix(0,pairsno,1)),
	as.double(matrix(0,pairsno,1)));
	score = d[[10]]/2;
	P = d[[11]];
	Q = d[[12]];
	retVal<-list(score=score,P=P, Q=Q,labels=labels,pairs=pairs);	

}



#KTSP.Train trains the KTSP classifier. 'data' is the matrix of the expression values 
# whose columns represents samples and the rows represents the genes. 
# 'situation' is a vector containing the training labels. Its elements should be one or zero.
# n is the number of top disjoint pairs. 
# If after picking some pairs, only pairs with score left,
# no more pair is chosen even the 'n' has not been reached. 

SWAP.KTSP.Train<-function( data, situation, krange= c(3,5:10), RestrictedPairs=c() , Filter=SWAP.Filter.Wilcoxon, FilterParam=list(UpDown=TRUE,geneNo=100),chooseKDyn=TRUE)
{
	
	if( chooseKDyn )
	{
		classifiers = SWAP.KTSP.Train( data, situation, krange, RestrictedPairs, Filter, FilterParam,chooseKDyn=FALSE);	
		classifierDynK = SWAP.ChoooseDynk( data, situation, classifiers, SWAP.KTSP.Statitsics);
		return( classifierDynK );

			
	}else{

		#Checks wether we should work on gene names or indices
		if(class(rownames(data))=="character")
			{
				usegenenames = TRUE;
			}else usegenenames = FALSE;
	
	
		
		#Filter gene
		if( length(Filter) > 0)
		{
			index = Filter(situation,data,FiltParam=FilterParam);
			if( (class(index) == "character" && usegenenames==FALSE) || (class(index) == "numeric" && usegenenames==TRUE) )
				stop("Filter's output must be consistent with gene names.") 
		}else{
			if( usegenenames )
			{
				index = rownames(data);
			}else
				index = 1:nrow(data);
		}

	
		if( length(RestrictedPairs)==0 )
		{
			classifiers = SWAP.KTSP.Train.Plain( data, situation, krange, index  );		
		}else{
			if( class(RestrictedPairs[1]) != class(index) )
				stop("RestictedPairs must be coherent with gene names."); 
			FilteredPairs = which(is.element(RestrictedPairs[,1],index)& is.element(RestrictedPairs[,2],index));
			scores = SWAP.CalculateSignedScore.Restricted(  situation , data, data, RestrictedPairs[FilteredPairs,]  );
			sortsc = sort(scores$score,decreasing=TRUE,index= TRUE);
			pairs = SWAP.KTSP.Pairs.Disjoint(RestrictedPairs[FilteredPairs [sortsc$ix],],
				sortsc$x,min(max(krange),length(scores)));
			maxTSPs = length(pairs$scores);
			classifiers = vector(mode="list",length=length(krange));

			for(i in 1:length(krange) )
			{
				classifiers[[i]]$name = sprintf('%dRTSP',krange[i]);
				classifiers[[i]]$TSPs = pairs$index[1:min(krange[i],maxTSPs),];
				classifiers[[i]]$score= pairs$scores[1:min(krange[i],maxTSPs),];
			}
		
		
		}
		return(classifiers);
	}
}

SWAP.KTSP.Train.Plain<-function( data, situation, krange= c(3,5:10), genes=c())
{
	
	
	situationV = as.vector(situation);
	if( length(unique(situationV)) != 2)
		stop("The situation must contain only exact two different variable.")

	if( length(genes)==0)
	{
		if( length(rownames(data))==0)
		{
			genes = 1:nrow(data);
		}else genes = rownames(data);
			
	}
	labels = unique(situationV);
	KTSPout =KTSP.CalculateSignedScore( situation , data[genes,], data[genes,] );
	
	TSPscore = vector(length=n);

	score = KTSPout$score;
	absscore = abs(score);
	maxK = max(krange);
	
	TSPsInd = matrix(0,maxK ,2);
	

	for(i in 1:maxK )
	{
		nextPair = arrayInd(which.max(absscore),.dim=dim(score));
		if( absscore[nextPair[1],nextPair[2]] < 1e-4 )
		{
			TSPsInd = TSPsInd[1:i,];
			TSPscore = TSPscore[1:i];
			break;
		}	

		if( score[nextPair[1],nextPair[2]]<0 )
		{
			TSPsInd[i,]= nextPair;
		}else
		 TSPsInd[i,]= c(nextPair[2],nextPair[1]);

		TSPscore[i] = absscore[TSPsInd[i,1],TSPsInd[i,2]]/2;
		
		absscore[TSPsInd[i,1],]=0;
		absscore[TSPsInd[i,2],]=0;
		absscore[,TSPsInd[i,1]]=0;
		absscore[,TSPsInd[i,2]]=0;

	}
		

	TSPs = cbind(genes[TSPsInd[,1]], genes[TSPsInd[,2]]);

	maxTSPs = nrow(TSPs);
	classifiers = vector(mode="list",length(krange));
	for(i in 1:length(krange) )
	{
		classifiers[[i]]$name = sprintf('%dTSP',krange[i]);
		classifiers[[i]]$TSPs = TSPs[1:min(krange[i],maxTSPs),];
		classifiers[[i]]$score = TSPscore[1:min(krange[i],maxTSPs)];		
		classifiers[[i]]$labels = labels; 		
	}
	return(classifiers);

}




SWAP.KTSP.Statitsics<-function( data, classifiers)
{
	if( length(names(classifiers))>0 )
	{
		if(length(classifiers$score)>2)		
		{
			KTSPstat =
			apply(data[classifiers$TSPs[,1],]>data[classifiers$TSPs[,2],],2,sum)-
			apply(data[classifiers$TSPs[,1],]<data[classifiers$TSPs[,2],],2,sum);		
		}else
		{		
			KTSPstat = apply(data[classifiers$TSPs[1],]>data[classifiers$TSPs[2],],2,sum)-
				     apply(data[classifiers$TSPs[1],]<data[classifiers$TSPs[2],],2,sum);			
		}

	}else
	{
		KTSPstat = matrix( 0, nrow = length(classifiers), ncol = ncol(data));

		for( i in 1:length(classifiers) )
		{

			if(length(classifiers[[i]]$score)>2)
			{				
				KTSPstat[i,] =
				apply(data[classifiers[[i]]$TSPs[,1],]>data[classifiers[[i]]$TSPs[,2],],2,sum)-
				apply(data[classifiers[[i]]$TSPs[,1],]<data[classifiers[[i]]$TSPs[,2],],2,sum);

			}else
			{
				KTSPstat[i,] =
				apply(data[classifiers[[i]]$TSPs[,1],]>data[classifiers[[i]]$TSPs[,2],],2,sum)-
				apply(data[classifiers[[i]]$TSPs[,1],]<data[classifiers[[i]]$TSPs[,2],],2,sum);

			}
		}
	}
	return(KTSPstat);
}

# The classifier for the test data. 'data' is the test data. 
SWAP.KTSP.Classify<-function( data, classifier,thresh=0)
{
	
	if( length(names(classifiers))>0 )
	{
		mylabels = classifier$labels;
	}else{
		mylabels = classifier[[1]]$labels;
	}
	ktspStat = SWAP.KTSP.Statitsics( data, classifier);
	pred = matrix( mylabels[1],nrow(ktspStat),ncol(ktspStat));
	classifiernames = vector(mode="character",length=length(classifier));
	
	
	for( i in 1:length(classifiers) )
	{
		pred[i,which(ktspStat[i,]>thresh)]= mylabels[1];
		pred[i,which(ktspStat[i,]<thresh)]= mylabels[2];
		classifiernames[i] = classifier[[i]]$name;
		
	}
	rownames(pred)<- classifiernames;
	return(pred);

}



# The classifier for the test data. 'data' is the test data. 

SWAP.Filter.Wilcoxon<-function(situation,data,FiltParam=list(genesNo=100,UpDown=TRUE))
{
	tiedData = SWAP.rank(data);
	tiedDataP = t(SWAP.rank(t(tiedData)));
	n = sum(situation==FALSE);
	m = sum(situation==TRUE);
	sumzeros = (apply(tiedDataP[,which(situation==situation[1])],1,sum))
	windex = (sumzeros -n*(n+m+1)/2)/sqrt(n*m*(n+m+1)/12);
	
		
	if( FiltParam$UpDown == TRUE )
	{	
		s = order(windex,decreasing=TRUE);
		lens = length(s);
		genesIndexUp = s[1:min(c(round(FiltParam$genesNo/2),lens))]; 
		genesIndexDown = s[ max(c(lens-round(FiltParam$genesNo/2),1)):lens]; 
		genesIndex =  unique(c(genesIndexUp,genesIndexDown));

	}else
	{
		s = order(abs(windex),decreasing=TRUE);
		genesIndex = s[1:min(c(FiltParam$genesNo,length(s)))]; 
	}

	if( length(rownames(data))==0)
	{
		return(genesIndex);
	}else 	
		return(rownames(data)[genesIndex]);
	
}

SWAP.ChoooseDynk<-function( data, situation, classifiers,Discriminant)
{
	s0 = which(situation==0);
	s1 = which(situation==1);


	
	mintt = -Inf;#minimum t-Test
	for( k in 1:length(classifiers) )
	{
		stat0 = Discriminant( data[,s0], classifiers[[k]]);
		stat1 = Discriminant( data[,s1], classifiers[[k]]);

		tt = abs(mean(stat0)-mean(stat1))/sqrt(var(stat1)+var(stat0)+0.000000001);#current t-test
		#If there is a tie score, we pick the classifier which shows up first
		if( abs(mintt-tt)>0.0000001 && mintt<tt)
		{
			mintt = tt;
			kmin = k;
		}		
		
	}
	return(classifiers[[kmin]]);
	
}

SWAP.KTSP.Pairs.Disjoint<-function(index,scores,k)
{
	sind = sort(scores,decreasing=TRUE,index=TRUE);
	pairssorted = index[sind$ix,];
	kTSPs = pairssorted[1,];
	scores = sind$x[1];
	n = 1;
	for( i in 2:nrow(index))
	{
		if(pairssorted[i,1]%in%kTSPs==FALSE &&pairssorted[i,2]%in%kTSPs==FALSE)
		{
			kTSPs = rbind(kTSPs,pairssorted[i,]);
			scores = rbind(scores,sind$x[i]);
			n = n + 1;
			if( n >= k )
				break;			
		}		
	}
		
	return(list(index=kTSPs,scores=scores));
}


