readLandmarks <- function(file, period.to.space = TRUE, row.names = 1, flip = FALSE, quote="\"", sep=',', header=TRUE){

	# Read file
	mat <- as.matrix(read.csv(file=file, row.names=row.names, quote=quote, sep=sep, header=header))

	# If name in columns, convert to array
	if(nchar(colnames(mat))[1] > 2){
		
		# Get names
		sp_names <- gsub('[_|.](|X|Y|Z)$', '', colnames(mat), ignore.case=TRUE)
		
		# Replace period with space
		if(period.to.space) sp_names <- gsub('[.]', ' ', sp_names)
		
		# Get unique names
		unique_sp_names <- unique(sp_names)

		# Detect whether 2D or 3D
		if(length(sp_names) == length(unique_sp_names)*3){ ndim <- 3 }else{ ndim <- 2 }

		# Create array
		if(flip){

			arr <- array(NA, dim=c(length(unique_sp_names), ndim, nrow(mat)), dimnames=list(unique_sp_names, c('x','y','z')[1:ndim], rownames(mat)))
		
			# Fill array
			for(unique_sp_name in unique_sp_names) arr[unique_sp_name, , ] <- t(mat[, unique_sp_name == sp_names])

		}else{

			arr <- array(NA, dim=c(nrow(mat), ndim, length(unique_sp_names)), dimnames=list(rownames(mat), c('x','y','z')[1:ndim], unique_sp_names))
		
			# Fill array
			for(unique_sp_name in unique_sp_names) arr[, , unique_sp_name] <- mat[, unique_sp_name == sp_names]
		}

		# Return array		
		return(arr)
	}

	mat
}