writeMotion <- function(x, file, digits = NULL, na="NaN", ...){

	# Create write mat
	write_mat <- NULL
	tmat_only <- FALSE

	if(is.matrix(x)){

		write_mat <- x

	}else{
	
		# If only object in list is tmat
		if(length(x) == 1 && names(x)[1] == 'tmat') tmat_only <- TRUE

		for(xn in names(x)){
	
			# Internal field used with filterRows
			if(xn %in% c('replace.rows', 'n.iter')) next

			if('tmat' %in% class(x[[xn]]) || xn == 'tmat'){
			
				# Convert array of transformation matrices to matrix
				tmat <- as.data.frame(tmarr2mat(x[[xn]]))

				if(!is.null(digits)) tmat <- signif(tmat, digits)

				# Add columns
				if(is.null(write_mat)){
					write_mat <- tmat
				}else{
					write_mat <- cbind(write_mat, tmat)
				}

			}else if('xyz' %in% class(x[[xn]]) || xn == 'xyz'){

				# Convert array of points to matrix
				xyz <- as.data.frame(arr2mat(x[[xn]]))
				
				if(!is.null(digits)) xyz <- signif(xyz, digits)

				# Add columns
				if(is.null(write_mat)){
					write_mat <- xyz
				}else{
					write_mat <- cbind(write_mat, xyz)
				}

			}else{
			
				# If is list, unlist
				if(is.list(x[[xn]])) x[[xn]] <- unlist(x[[xn]])

				# Need to test whether values are numeric first (including case of some or all NAs)
				if(!is.null(digits) && sum(!suppressWarnings(is.numeric(x[[xn]]))) == 0 && sum(is.na(x[[xn]])) == 0) x[[xn]] <- signif(x[[xn]], digits)

				# Add columns
				if(is.null(write_mat)){
					write_mat <- as.data.frame(matrix(x[[xn]], ncol=1, dimnames=list(NULL, xn)))
				}else{
					write_mat <- cbind(write_mat, as.data.frame(matrix(x[[xn]], ncol=1, dimnames=list(NULL, xn))))
				}
			}
		}
	}
	
	# If only tmat set NAs to NaN so that the file is read properly by Maya?
	#print(tmat_only)
	#if(tmat_only){ na_handling <- "NaN" }else{ na_handling <- "NA" }

	if(grepl('[.]csv$', file)) write.csv(x=write_mat, file=file, row.names=FALSE, na=na, ...)
}