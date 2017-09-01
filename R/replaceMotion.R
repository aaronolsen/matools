replaceMotion <- function(motion, replacement, replace.rows = NULL){

	# Check for
	if(is.null(replace.rows) && !'replace.rows' %in% names(replacement)) stop("'replace.rows' not specified in function call or found in 'replacement' object. Rows to be replaced must be specified.")

	# Get number of rows
	if(!is.null(motion$n.iter)){

		num_rows <- motion$n.iter

	}else{

		for(xn in names(motion)){

			if(xn %in% c('replace.rows', 'n.iter')) next

			if(is.null(dim(motion[[xn]]))){
				num_rows <- length(motion[[xn]])
			}else{
				if(length(dim(motion[[xn]])) == 2){
					num_rows <- dim(motion[[xn]])[2]
				}else if(length(dim(motion[[xn]])) == 3){
					num_rows <- dim(motion[[xn]])[3]
				}else if(length(dim(motion[[xn]])) == 4){
					num_rows <- dim(motion[[xn]])[4]
				}
			}
		}
	}
	
	for(xn in names(replacement)){

		# Skip if replace.rows
		if(xn %in% c('replace.rows', 'n.iter')) next

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
			motion[[xn]][replacement$replace.rows] <- replacement[[xn]]
		}else{
			if(length(dim_rows) == 2){
				motion[[xn]][, replacement$replace.rows] <- replacement[[xn]]
			}else if(length(dim_rows) == 3){
				motion[[xn]][, , replacement$replace.rows] <- replacement[[xn]]
			}else if(length(dim_rows) == 4){
				motion[[xn]][, , , replacement$replace.rows] <- replacement[[xn]]
			}
		}
	}

	motion
}
