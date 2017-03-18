readTransform <- function(file, source='xmalab'){

	# Check if file or directory
	if(grepl('[.]csv$', file, ignore.case=TRUE)){

	}else{

		# Add if forward slash to end if not present
		if(!grepl('/$', file)) file <- paste0(file, '/')

		# Get file list
		list_files <- list.files(file)
		
		# Get names of rigid bodies
		body_names <- gsub('(RigidBody[0-9]+_)|_transformation[.]csv$', '', list_files)
		
		# Get number of frames from first file
		nframes <- dim(as.matrix(read.csv(paste0(file, list_files[1]))))[1]
		
		# Create array for all transformations
		tr_arr <- array(NA, dim=c(4,4,nframes,length(body_names)), dimnames=list(NULL, NULL, NULL, body_names))
		
		# Fill transformation array
		for(i in 1:length(list_files)){

			# Read csv
			read_csv <- as.matrix(read.csv(paste0(file, list_files[i])))
			
			# Add to array
			tr_arr[, , , i] <- t(read_csv)
		}
	}

	tr_arr
}
