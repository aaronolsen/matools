removeRows <- function(motion, remove.rows){

	# Get number of rows
	num_rows <- nrowMotion(motion)

	# Set which rows to keep
	select_which <- rep(TRUE, num_rows)
	select_which[remove.rows] <- FALSE
	
	# Get indices of removed rows from original motion object if available
	if(!is.null(motion$replace.rows)) replace_rows <- motion$replace.rows[remove.rows]

	for(motion_name in names(motion)){

		# Internal field used with subsetRows
		if(motion_name %in% c('n.iter')) next
		
		motion_class <- class(motion[[motion_name]])

		dim_motion <- dim(motion[[motion_name]])
		
		if(is.null(dim_motion)){
			motion[[motion_name]] <- motion[[motion_name]][select_which]
		}else{
			if(length(dim_motion) == 2){
				#if(dim_motion[1] == initial_length) motion[[motion_name]] <- motion[[motion_name]][select_which, ]
				motion[[motion_name]] <- motion[[motion_name]][, select_which]
			}else if(length(dim_motion) == 3){
				motion[[motion_name]] <- motion[[motion_name]][, , select_which]
			}else if(length(dim_motion) == 4){
				motion[[motion_name]] <- motion[[motion_name]][, , , select_which]
			}
		}

		class(motion[[motion_name]]) <- motion_class
	}

	# Save where rows were removed in original motion object in case they will be replaced
	# remove.rows only needed if motion object is a subset (replace.rows is not NULL)
	if(!is.null(motion$replace.rows)){
		
		# Add removed rows to remove.rows, using original motion object indices
		if(is.null(motion$remove.rows)){
			motion$remove.rows <- replace_rows
		}else{
			motion$remove.rows <- c(motion$remove.rows, replace_rows)
		}
	}
	
	# Set new number of rows/iterations
	motion$n.iter <- sum(select_which == TRUE)
	
	motion
}
