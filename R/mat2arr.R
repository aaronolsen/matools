mat2arr <- function(mat, pattern = '[.|_](x|y|z)$', ignore.case = TRUE){

	if(pattern == '[.|_](x|y)$'){

		# Create array
		arr <- array(NA, dim=c(ncol(mat)/2,2,nrow(mat)), 
			dimnames=list(unique(gsub(pattern, '', colnames(mat), ignore.case=ignore.case)), NULL, rownames(mat)))

		# Question marks at end of marker names become '.'
		dimnames(arr)[[1]] <- gsub('[.]$', '', dimnames(arr)[[1]])

		# Fill array
		for(i in 1:dim(arr)[3]) arr[, , i] <- matrix(mat[i, ], nrow=ncol(mat)/2, ncol=2, byrow=TRUE)

	}else{

		# Create array
		arr <- array(NA, dim=c(ncol(mat)/3,3,nrow(mat)), 
			dimnames=list(unique(gsub(pattern, '', colnames(mat), ignore.case=ignore.case)), NULL, rownames(mat)))

		# Question marks at end of marker names become '.'
		dimnames(arr)[[1]] <- gsub('[.]$', '', dimnames(arr)[[1]])

		# Fill array
		for(i in 1:dim(arr)[3]) arr[, , i] <- matrix(mat[i, ], nrow=ncol(mat)/3, ncol=3, byrow=TRUE)
	}
	
	arr
}
