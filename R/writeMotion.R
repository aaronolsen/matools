writeMotion <- function(x, file){

	# Create write mat
	write_mat <- NULL

	for(xn in names(x)){
	
		if('tmat' %in% class(x[[xn]])){
			
			# Convert array of transformation matrices to matrix
			tmat <- tmarr2mat(x[[xn]])
			
			# Add columns
			write_mat <- cbind(write_mat, tmat)

		}else if('xyz' %in% class(x[[xn]])){

			# Convert array of points to matrix
			xyz <- arr2mat(x[[xn]])
			
			# Add columns
			write_mat <- cbind(write_mat, xyz)

		}else{

			# Add columns
			write_mat <- cbind(write_mat, matrix(x[[xn]], ncol=1, dimnames=list(NULL, xn)))
		}
	}
	
	if(grepl('[.]csv$', file)) write.csv(x=write_mat, file=file)
	
	NULL
}