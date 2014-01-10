##################################################
##################################################
### This is the plain KTSP training function with only genes
### Using all the feature provided in data (i.e. all the rows)
### Filtering can be achieved by subsetting data BEFORE training
### Example 1: SWAP.KTSP.Train.Plain(matETABM115, ttmGroupETABM115)


### Situation perhaps should be a factor! Make sure to change this

SWAP.KTSP.Train.Plain <- function(data, situation, krange= c(3:10)) {

	## GET THE GENES
	genes <- rownames(data)

	## Verify outcome lables
	situationV <- as.vector(situation);

	## Define the maximun number of TSP
	maxK <- max(krange);

	## Evaluate the argument 'situation'
	if (length(unique(situationV)) != 2) {
		stop("The situation must contain only exact two different variable.")
	}

	## Get the class lables
	labels <- unique(situationV);

	## Compute the scores
	KTSPout <- SWAP.CalculateSignedScore(situation , data, data)

	## Process the score
	TSPscore <- vector(length=maxK)
	score <- KTSPout$score
	absscore <- abs(score)

	## Get the index for the scores
	TSPsInd <- matrix(0, maxK, 2)

	for (i in 1:maxK) {
		nextPair <- arrayInd(which.max(absscore), .dim=dim(score))

		## If the score of the pairs are negligible do not go down the list
		if (absscore[ nextPair[1], nextPair[2] ] < 1e-4 ) {
			TSPsInd <- TSPsInd[ 1:i, ]
			TSPscore <- TSPscore[ 1:i]
			break
		}

		if (score[nextPair[1], nextPair[2]]<0) {
			TSPsInd[i,]= nextPair
		} else {
			TSPsInd[i,] = c(nextPair[2], nextPair[1])
		}

		TSPscore[i] <- absscore[TSPsInd[i,1], TSPsInd[i,2]]/2
		absscore[ TSPsInd[i , 1] , ] <- 0
		absscore[ TSPsInd[i , 2] , ] <- 0
		absscore[ , TSPsInd[i , 1] ] <- 0
		absscore[ , TSPsInd[i , 2] ] <- 0
	}

	## PREPARE TSP
	TSPs <- cbind(genes[TSPsInd[ , 1]]  , genes[TSPsInd[ , 2]])

	maxTSPs <- nrow(TSPs)
	classifiers <- vector(mode="list", length(krange))

	## SELECT THE KTSP IN THE K-RANGE
	for (i in 1:length(krange)) {
		classifiers[[i]]$name = sprintf('%dTSP',krange[i])

		if (min(krange[i], maxTSPs) == 1) {
			classifiers[[i]]$TSPs <- matrix( TSPs[1,], ncol=2, nrow=1)
		} else {
			classifiers[[i]]$TSPs <- TSPs[1:min(krange[i], maxTSPs),]
		}
		classifiers[[i]]$score <- TSPscore[1:min(krange[i], maxTSPs)]
		classifiers[[i]]$labels <- labels
	}

	## Reeturn the classifiers
	return(classifiers)

}


##################################################
##################################################
### SWAPKTSP.Train trains the KTSP classifier.
### 'data' is the matrix of the expression values
### whose columns represents samples and the rows represents the genes.
### ROWNAMES WILL BE USED AS GENE NAMES
### 'situation' is a vector containing the training labels.
### 'RestrictedPairs' is a list: one vector of length 2 for each TSP

### IS THE FOLLOWING STILL TRUE?
### Its elements should be one or zero.
### n is the number of top disjoint pairs.
### If after picking some pairs, only pairs with score left,
### no more pair is chosen even the 'n' has not been reached.

### Bahman, I think, and modified accordingly the code:
### 1) We should force to use rownames as gene names;
### 2) No FILTER! People can do this upfront on the matrix, before starting training.

SWAP.KTSP.Train <- function(data, situation,
			    krange = c(1:10),
			    RestrictedPairs=NULL,
			    chooseKDyn = TRUE) {

	## Create the index using rownames(data)
	index <- rownames(data)

	## Test for the availability of restricted pairs
	if (is.null(RestrictedPairs)) {
		classifiers <- SWAP.KTSP.Train.Plain(data, situation, krange, index)
	} else {
		## Test restricted gene pairs availability in data
			## First check if both features in a TSP are available
			keepPairs <- sapply(RestritedPairs, function(x) any(! x %in% rownames(data)) )
			RestrictedPairs <- RestrictedPairs[keepPairs]

			## Stop is none of the TSP survived filtering
			if (length(RestrictedPairs) == =) {
				stop("None of the 'RestictedPairs' is available in 'data'");
			}

			## Compute and sort the scores for the restricted pairs
			scores <- SWAP.CalculateSignedScore.Restricted(situation, data, data,
								       RestrictedPairs)
			sortsc <- sort(scores$score, decreasing=TRUE, index= TRUE)

			## Compute disjoint TSPs
			pairs <- SWAP.KTSP.Pairs.Disjoint(RestrictedPairs[sortsc$ix],
							  sortsc$x,min(max(krange), length(scores)))
			maxTSPs <- length(pairs$scores)
			classifiers <- vector(mode="list", length=length(krange))

			## Prepare the output for each TSP
			for (i in 1:length(krange)) {
				classifiers[[i]]$name <- sprintf('%dRTSP', krange[i]);
				classifiers[[i]]$TSPs <- pairs$index[1:min(krange[i], maxTSPs),]
				classifiers[[i]]$score <- pairs$scores[1:min(krange[i], maxTSPs),];
			}
		}

	## Chose the best K TSP
	if (chooseKDyn) {
		classifiers <- SWAP.ChoooseDynk(data, situation, classifiers, SWAP.KTSP.Statitsics)
	}

	## Return the classifier
	return(classifiers)
}



##################################################
##################################################
### It calculates the KTSP statistics : \sum_{k} I(X_i_k <X_j_k) - I(X_j_k <X_i_k)
### The usage is:
### statistics <- SWAP.KTSP.Statitsics(matETABM115, classifiers)

SWAP.KTSP.Statitsics <- function(data, classifiers) {

	## Are there classifier available?
	if (length(names(classifiers))>0) {
		## Is 'data' a vector?
		if (is.vector(data)) {
			KTSPstat <- sum(data[ classifiers$TSPs[,1]] > data[classifiers$TSPs[,2]])
			-
				sum(data[classifiers$TSPs[,1]] < data[classifiers$TSPs[,2]])
			KTSPstat <- matrix(KTSPstat, nrow = 1, ncol=1);
		} else {
			## If  'data' is not a vector
			KTSPstat <- apply(data[classifiers$TSPs[,1],] > data[classifiers$TSPs[,2],], 2, sum)
			-
				apply(data[classifiers$TSPs[,1],] < data[classifiers$TSPs[,2],],2,sum)
			KTSPstat <- matrix(KTSPstat,nrow = 1);
		}

		##Name of the classifier
		rownames(KTSPstat) <- classifiers$name

	## If there are not  available classifiers do the following?
	} else {
		## Is 'data' a vector?
		if (is.vector(data)) {
			KTSPstat <- matrix( 0, nrow = length(classifiers), ncol = 1);
		} else {
			## If  'data' is not a vector
			KTSPstat <- matrix( 0, nrow = length(classifiers), ncol = ncol(data))
		}

		## Classifier names
		classifiernames <- vector(mode="character", length=length(classifiers))

		## For each classifier fo the following
		for (i in 1:length(classifiers))	{
			## If 'data' is a vector
			if (is.vector(data)) {
				KTSPstat[i] <-	sum(data[classifiers[[i]]$TSPs[,1]] > data[classifiers[[i]]$TSPs[,2]])
				-
					sum(data[classifiers[[i]]$TSPs[,1]] < data[classifiers[[i]]$TSPs[,2]])
			} else {
				## If 'data' is not a vector
				KTSPstat[i,] <- apply(data[classifiers[[i]]$TSPs[,1],] > data[classifiers[[i]]$TSPs[,2],], 2 ,sum)
				-
					apply(data[classifiers[[i]]$TSPs[,1],] < data[classifiers[[i]]$TSPs[,2],], 2, sum)
			}

			## Classifiers names
			classifiernames[i] <- classifiers[[i]]$name;
		}

		## Assign the classifier names
		rownames(KTSPstat) <- classifiernames;
	}

	### Return the results
	return(KTSPstat)
}


##################################################
##################################################
### The classifier for the test data. 'data' is the test data.
### just threshold the outcome of KTSP.Statistics

### The argument 'thres' is NOT USED, I will take it aout
### SWAP.KTSP.Classify <- function(data, classifiers, thresh=0) {

SWAP.KTSP.Classify <- function(data, classifiers) {

	##Checks if there is only one kTSP or more
	if (length(names(classifiers)) > 0) {
		mylabels <- classifiers$labels
	}

	## Classify with the data
	ktspStat <- SWAP.KTSP.Statitsics(data, classifiers)

	## Compute tha return the classification result
	out <- ifelse(ktspStat > 0, mylabels[[2]], mylabels[[1]])
	return(out)
}


##################################################
##################################################
### We may need to write this function in C
### Not yet exported in the NAMESPACE

SWAP.KTSP.Pairs.Disjoint <- function(index, scores, k) {

	sind <- sort(scores, decreasing=TRUE, index=TRUE)
	pairssorted <- index[sind$ix,]

	kTSPs <- pairssorted[1,]
	scores <- sind$x[1]

	n <- 1

	for (i in 2:nrow(index)) {
		if (pairssorted[i,1] %in% kTSPs == FALSE && pairssorted[i,2] %in% kTSPs==FALSE) {
			kTSPs <- rbind(kTSPs, pairssorted[i,])
			scores <- rbind(scores,sind$x[i])
			n <- n + 1

			if (n >= k) {
				break
			}
		}
	}

	out <- list(index=kTSPs, scores=scores)
	return(out)

}


