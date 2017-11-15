writeMotion <- function(x, file){

	# Create write mat
	write_mat <- NULL

	if(is.matrix(x)){

		write_mat <- x
	}else{

		for(xn in names(x)){
	
			# Internal field used with filterRows
			if(xn %in% c('replace.rows', 'n.iter')) next

			if('tmat' %in% class(x[[xn]]) || xn == 'tmat'){
			
				# Convert array of transformation matrices to matrix
				tmat <- as.data.frame(tmarr2mat(x[[xn]]))
			
				# Add columns
				write_mat <- cbind(write_mat, tmat)

				# Add columns
				if(is.null(write_mat)){
					write_mat <- tmat
				}else{
					write_mat <- cbind(write_mat, tmat)
				}

			}else if('xyz' %in% class(x[[xn]]) || xn == 'xyz'){

				# Convert array of points to matrix
				xyz <- as.data.frame(arr2mat(x[[xn]]))

				# Add columns
				if(is.null(write_mat)){
					write_mat <- xyz
				}else{
					write_mat <- cbind(write_mat, xyz)
				}

			}else{
			
				# If is list, unlist
				if(is.list(x[[xn]])) x[[xn]] <- unlist(x[[xn]])

				# Add columns
				if(is.null(write_mat)){
					write_mat <- as.data.frame(matrix(x[[xn]], ncol=1, dimnames=list(NULL, xn)))
				}else{
					write_mat <- cbind(write_mat, as.data.frame(matrix(x[[xn]], ncol=1, dimnames=list(NULL, xn))))
				}
			}
		}
	}

	if(grepl('[.]csv$', file)) write.csv(x=write_mat, file=file, row.names=FALSE)
	
	NULL
}