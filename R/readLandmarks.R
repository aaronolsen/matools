readLandmarks <- function(file, period.to.space = TRUE, row.names = 1, flip = FALSE, quote="\"", sep=',', header=TRUE){

	if(!file.exists(file)) stop(paste0("File '", file, "' not found."))

	# Read file
	mat <- suppressWarnings(as.matrix(read.csv(file=file, row.names=row.names, quote=quote, sep=sep, header=header)))

	# Single row of points, reformat as matrix
	if(dim(mat)[1] == 1){
		
		# Read file
		mat <- suppressWarnings(as.matrix(read.csv(file=file, row.names=NULL, quote=quote, sep=sep, header=header)))

		# Get landmark names
		sp_names <- unique(gsub('[_|.](|x|y|z)$', '', colnames(mat), ignore.case=TRUE))
		
		# Create matrix
		as_mat <- matrix(mat, nrow=dim(mat)[2]/3, ncol=3, dimnames=list(sp_names, NULL), byrow=TRUE)
		
		return(as_mat)
	}

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