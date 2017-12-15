subsetRows <- function(motion, criteria, new = FALSE){

	for(criteria_name in names(criteria)){
	
		if(!criteria_name %in% names(motion)) stop("'", criteria_name, "' not found in motion property names.")
	
		#
		if(length(criteria[[criteria_name]]) > 1){
			select_which <- which(motion[[criteria_name]] %in% criteria[[criteria_name]])
		}else{
			select_which <- which(grepl(criteria[[criteria_name]], motion[[criteria_name]]))
		}
		
		#
		#initial_length <- length(criteria[[criteria_name]])
		
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

		# Save where rows were extracted in case they will be replaced
		if(is.null(motion$replace.rows)){
			motion$replace.rows <- select_which
		}else{
			motion$replace.rows <- motion$replace.rows[select_which]
		}
	}
	
	# Set new number of rows/iterations
	motion$n.iter <- length(motion$replace.rows)

	if(new) motion$replace.rows <- NULL
	
	motion
}