filterRows <- function(rows, criteria){

	for(criteria_name in names(criteria)){
	
		if(!criteria_name %in% names(rows)) stop("'", criteria_name, "' not found in rows names.")
	
		#
		select_which <- which(rows[[criteria_name]] %in% criteria[[criteria_name]])
		
		#
		initial_length <- length(criteria[[criteria_name]])
		
		for(rows_name in names(rows)){
			
			rows_class <- class(rows[[rows_name]])

			dim_rows <- dim(rows[[rows_name]])
			
			if(is.null(dim_rows)){
				rows[[rows_name]] <- rows[[rows_name]][select_which]
			}else{
				
				if(length(dim_rows) == 2){
					if(dim_rows[1] == initial_length) rows[[rows_name]] <- rows[[rows_name]][select_which, ]
					if(dim_rows[2] == initial_length) rows[[rows_name]] <- rows[[rows_name]][, select_which]
				}else if(length(dim_rows) == 3){
					rows[[rows_name]] <- rows[[rows_name]][, , select_which]
				}else if(length(dim_rows) == 4){
					rows[[rows_name]] <- rows[[rows_name]][, , , select_which]
				}
			}

			class(rows[[rows_name]]) <- rows_class
		}
	}

	rows
}
