replaceRows <- function(motion, replacement, replace.rows = NULL, remove.rows = NULL){

	# Check for
	if(is.null(replace.rows) && !'replace.rows' %in% names(replacement)) stop("'replace.rows' not specified in function call or found in 'replacement' object. Rows to be replaced in 'motion' must be specified.")

	#
	if(!is.null(motion$replace.rows)) warning("'replace.rows' is present in 'motion' which means you are replacing rows in an object that is a subset with rows that are also a subset. Rows may not be replaced in the correct order.")

	if(!is.null(replace.rows)){
		replace_rows <- replace.rows
	}else{
		replace_rows <- replacement[['replace.rows']]
	}

	remove_rows <- NULL
	if(!is.null(remove.rows)){
		remove_rows <- remove.rows
	}else{
		if(!is.null(replacement[['remove.rows']])) remove_rows <- replacement[['remove.rows']]
	}

	# Get number of rows
	num_rows <- nrowMotion(motion)
	
	# Set which rows to keep (rows only removed if remove_rows is non NULL)
	keep_rows <- rep(TRUE, num_rows)
	if(!is.null(remove_rows)) keep_rows[remove_rows] <- FALSE
	
	for(xn in names(replacement)){

		# Skip if replace.rows
		if(xn %in% c('replace.rows', 'remove.rows', 'n.iter')) next

		# If new motion has object not in original, add new object
		if(!xn %in% names(motion)){

			if(xn == 'tmat'){
				motion[[xn]] <- array(NA, dim=c(dim(replacement[[xn]])[1:3], num_rows), dimnames=list(NULL, NULL, dimnames(replacement[[xn]])[[3]], NULL))
			}else if(xn == 'xyz'){
				motion[[xn]] <- array(NA, dim=c(dim(replacement[[xn]])[1:2], num_rows), dimnames=list(dimnames(replacement[[xn]])[[1]], NULL, NULL))
			}else{
				motion[[xn]] <- rep(NA, length=num_rows)
			}
		}
	
		dim_rows <- dim(motion[[xn]])
		
		if(is.null(dim_rows)){
			
			motion[[xn]][replace_rows] <- replacement[[xn]]
			if(!is.null(remove_rows)) motion[[xn]] <- motion[[xn]][keep_rows]
			
		}else{
			if(length(dim_rows) == 2){
				motion[[xn]][, replace_rows] <- replacement[[xn]]
				if(!is.null(remove_rows)) motion[[xn]] <- motion[[xn]][, keep_rows]
			}else if(length(dim_rows) == 3){
			
				if(!is.null(dimnames(motion[[xn]]))){

					# Names in replacement not found in motion
					if(sum(!dimnames(replacement[[xn]])[[1]] %in% dimnames(motion[[xn]][, , replace_rows])[[1]]) > 0){
						
						# Both names
						both_names <- unique(c(dimnames(replacement[[xn]])[[1]], dimnames(motion[[xn]][, , replace_rows])[[1]]))
						
						# Create new array
						new_array <- array(NA, dim=c(length(both_names), dim(motion[[xn]])[2:3]), dimnames=list(both_names, 
							dimnames(motion[[xn]])[[2]], dimnames(motion[[xn]])[[3]]))
						
						# Fill new array
						new_array[dimnames(motion[[xn]])[[1]], , ] <- motion[[xn]]
						motion[[xn]] <- new_array
					}
					
					motion[[xn]][dimnames(replacement[[xn]])[[1]], , replace_rows] <- replacement[[xn]]

				}else{
					motion[[xn]][, , replace_rows] <- replacement[[xn]]
				}
			
				if(!is.null(remove_rows)) motion[[xn]] <- motion[[xn]][, , keep_rows]
			}else if(length(dim_rows) == 4){
				motion[[xn]][, , , replace_rows] <- replacement[[xn]]
				if(!is.null(remove_rows)) motion[[xn]] <- motion[[xn]][, , , keep_rows]
			}
		}
	}
	
	motion$n.iter <- sum(keep_rows)

	motion
}