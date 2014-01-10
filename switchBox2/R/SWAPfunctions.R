##################################################
##################################################
### Not exported in the NAMESPACE

SWAP.rank <- function(x) {
	return((apply(x,2,rank)));
}



##################################################
##################################################

SWAP.CalculateSignedScore <- function(situation , data1, data2) {

	n <- length(situation);
	m1 <- nrow(data1);
	m2 <- nrow(data2);

	datatied <- SWAP.rank(rbind(data1,data2));
	data1tied <- datatied[1:m1,];
	data2tied <- datatied[(m1+1):(m1+m2),];

	situationV <- as.vector(situation);

	if (length(unique(situationV)) != 2) {
		stop("The situation must contain only exact two different variable.")
	}

	labels <- unique(situationV);

	d <-	.C(
		"CalculateSignedScoreCore",
		as.integer(as.numeric(situationV==labels[1])), as.integer(n),
		as.double(data1tied), as.integer(m1),
		as.double(data2tied), as.integer(m2),
		as.double(matrix(0, m1, m2)),
		as.double(matrix(0, m1, m2)),
		as.double(1),
		as.double(matrix(0, m1, m2)),
		as.double(1),
		as.double(matrix(0, m1, m2)),
		as.double(matrix(0, m1, m2))
		);

	score <- (matrix(d[[7]], nrow=m1))/2;
#	a <- (matrix(d[[7]], nrow=m1));
#     M <- (matrix(d[[8]], nrow=1));
#	k <- (matrix(d[[9]], nrow=m1));
#	N <- (matrix(d[[10]], nrow=1));
     	P <- (matrix(d[[12]], nrow=m1));
	Q <- (matrix(d[[13]], nrow=m1));

	if (length(rownames(data1))>0 &&  length(rownames(data2))>0 ) {
		names1 <- rownames(data1);
		names2 <- rownames(data2);

		rownames(score) <- names1;
		colnames(score) <- names2;

		rownames(P) <- names1;
		colnames(P) <- names2;

		rownames(Q) <- names1;
		colnames(Q) <- names2;
	}

	retVal <- list(score=score, P=P, Q=Q, labels=labels);

}


##################################################
##################################################

SWAP.CalculateSignedScore.Restricted <- function(situation , data1, data2, pairs) {

	n <- length(situation);
	m1 <- nrow(data1);
	m2 <- nrow(data2);
	pairsno <- nrow(pairs);

	datatied <- SWAP.rank(rbind(data1, data2));
	data1tied <- datatied[1:m1,];
	data2tied <- datatied[(m1+1):(m1+m2),];

	situationV <- as.vector(situation);

	if (length(unique(situationV)) != 2) {
		stop("The situation must contain only exact two different variable.")
	}

	labels <- unique(situationV);


	if (class(pairs[1,1])=="character") {
		pairsind <- matrix(0,pairsno ,2);
		pairsind[,1] <- match(pairs[,1],rownames(data1));
		pairsind[,2] <- match(pairs[,2],rownames(data2));

		if (sum(is.na(pairsind))>0) {
			stop("The genes in pairs do not show in data1 or data2.");
		}
	} else {
		if (min(pairs)<1 || max(pairs[,1])>nrow(data1) || max(pairs[,2])>nrow(data2)) {
			stop("The genes in pairs do not show in data1 or data2.");
		}
		pairsind <- pairs;
	}

	d <-	.C(
		"CalculateSignedScoreRestrictedPairsCore",
		as.integer(as.numeric(situationV==labels[1])),
		as.integer(n),
		as.double(data1tied), as.integer(m1),
		as.double(data2tied), as.integer(m2),
		as.integer(pairsind[,1]-1),
		as.integer(pairsind [,2]-1),
		as.integer(pairsno),
		as.double(matrix(0, pairsno, 1)),
		as.double(matrix(0, pairsno, 1)),
		as.double(matrix(0, pairsno, 1))
		);

	score <- d[[10]]/2;
	P <- d[[11]];
	Q <- d[[12]];

	retVal <- list(score=score, P=P, Q=Q, labels=labels, pairs=pairs);

}




##################################################
##################################################
### The classifier for the test data. 'data' is the test data.

SWAP.Filter.Wilcoxon <- function(situation, data,
				 FiltParam=list(genesNo=100, UpDown=TRUE) ) {

	tiedData <- SWAP.rank(data);
	tiedDataP <- t(SWAP.rank(t(tiedData)));
	n <- sum(situation==FALSE);
	m <- sum(situation==TRUE);
	sumzeros <- (apply(tiedDataP[,which(situation==situation[1])],1,sum))
	windex <- (sumzeros -n*(n+m+1)/2)/sqrt(n*m*(n+m+1)/12);

	if (FiltParam$UpDown == TRUE) {
		s <- order(windex,decreasing=TRUE);
		lens <- length(s);
		genesIndexUp <- s[1:min(c(round(FiltParam$genesNo/2),lens))];
		genesIndexDown <- s[ max(c(lens-round(FiltParam$genesNo/2),1)):lens];
		genesIndex <- unique(c(genesIndexUp,genesIndexDown));
	} else {
		s <- order(abs(windex),decreasing=TRUE);
		genesIndex <- s[1:min(c(FiltParam$genesNo,length(s)))];
	}

	if (length(rownames(data))==0) {
		return(genesIndex);
	} else {
		return(rownames(data)[genesIndex]);
	}
}


##################################################
##################################################

SWAP.ChoooseDynk <- function(data, situation, classifiers, Discriminant) {

	s0 <- which(situation==0);
	s1 <- which(situation==1);

	##if there is only one classifier, return that classifier
	if (length(names(classifiers))>0) {
		return(classifiers);
	}

	##minimum t-Test
	mintt <- -Inf;
	for (k in 1:length(classifiers)) {
		stat0 <- as.vector(Discriminant( data[,s0], classifiers[[k]]));
		stat1 <- as.vector(Discriminant( data[,s1], classifiers[[k]]));

		##current t-test
		tt <- abs(mean(stat0)-mean(stat1))/sqrt(var(stat1)+var(stat0)+0.000000001);

		##If there is a tie score, we pick the classifier which shows up first
		if (abs(mintt-tt)>0.0000001 && mintt<tt) {
			mintt = tt;
			kmin = k;
		}
	}
	return(classifiers[[kmin]]);

}



