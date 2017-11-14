arr2mat <- function(arr, xyz=c('_X', '_Y', '_Z'), dim.to.row=3){

	# Get dimensions
	dim_arr <- dim(arr)

	# Dimnames
	dimnames_arr <- dimnames(arr)

	#print(paste0(lm_repeat, '.', xyz))
	#print(list(dimnames_arr[[3]], paste0(dimnames_arr[[1]], '.', letters[24:26])))
	
	# Convert to matrix
	if(dim.to.row == 3){

		# Repeate first dimension names
		lm_repeat <- c(t(matrix(dimnames_arr[[1]], nrow=length(dimnames_arr[[1]]), ncol=3)))

		# Create matrix
		mat <- matrix(NA, nrow=dim_arr[3], ncol=dim_arr[1]*dim_arr[2], 
			dimnames=list(NULL, paste0(lm_repeat, xyz)))
	
		# Fill matrix
		for(i in 1:nrow(mat)) mat[i, ] <- t(arr[, , i])

	}else if(dim.to.row == 1){

		# Repeate third dimension names
		lm_repeat <- c(t(matrix(dimnames_arr[[3]], nrow=length(dimnames_arr[[3]]), ncol=3)))

		# Create matrix
		mat <- matrix(NA, nrow=dim_arr[1], ncol=dim_arr[2]*dim_arr[3], 
			dimnames=list(dimnames_arr[[1]], paste0(lm_repeat, xyz)))
	
		# Fill matrix
		for(i in 1:nrow(mat)) mat[i, ] <- arr[i, , ]
	}

	mat
}
