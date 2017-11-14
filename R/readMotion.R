readMotion <- function(file, nrows = -1){

	## Reads in matrix of coordinates over time, with or without time column, or 
	## transformation matrices. File type is detected based on whether first column name 
	## ends in R11

	# Read first two lines
	read_lines <- readLines(file, n=2)

	# Check for quotes
	num_quotes <- length(gregexpr('"', read_lines[1])[[1]])
	num_commas <- c(length(gregexpr(',', read_lines[1])[[1]]), length(gregexpr(',', read_lines[2])[[1]]))

	# Check for ' , ' separators
	if(grepl(' , ', read_lines[2])){ sep <- ' , ' }else{ sep <- ',' }

	# Split first line at commas
	comma_split <- strsplit(read_lines[1], ',')[[1]]

	# Set row.names parameter
	#if(comma_split[1] == '\"\"' || num_commas[1] != num_commas[2]){ row.names <- 1 }else{ row.names <- NULL }

	# Set quote parameter
	if(num_quotes == num_commas[1]*2 + 2){ quote <- '\"' }else{ quote <- '' }

	# Read file
	if(grepl('[.]csv$', file[1])){
		read_matrix <- as.matrix(read.csv(file[1], nrows=nrows))
	}else if(grepl('[.]txt$', file[1])){
		stop("'.txt' file reading not yet supported.")
	}

	# If first column is X, remove
	row_names <- NULL
	if(colnames(read_matrix)[1] == 'X'){
		row_names <- read_matrix[, 1]
		read_matrix <- read_matrix[, 2:ncol(read_matrix)]
	}

	# Sort column names
	read_matrix <- read_matrix[, sort(colnames(read_matrix))]

	# Check if there are transformations
	tmat_cols <- grepl('_(R[1-3]{2}|[0-3]{2}|1|TX|TY|TZ)$', colnames(read_matrix), ignore.case=TRUE)

	# Check if there are xyz coordinates
	xyz_cols <- grepl('[_|.](|X|Y|Z)$', colnames(read_matrix), ignore.case=TRUE)

	# Check for additional info columns
	info_cols <- tmat_cols+xyz_cols == 0

	# Set logicals
	has_tmat <- ifelse(sum(tmat_cols) > 0, TRUE, FALSE)
	has_xyz <- ifelse(sum(xyz_cols) > 0, TRUE, FALSE)
	has_info <- ifelse(sum(info_cols) > 0, TRUE, FALSE)

	# Fill transformation array
	if(has_tmat){

		# Get tm matrix
		tm_matrix <- read_matrix[, tmat_cols]

		# Transformation matrix
		tmat_mat <- matrix(suppressWarnings(as.numeric(tm_matrix)), nrow(tm_matrix), ncol(tm_matrix), 
			dimnames=dimnames(tm_matrix))
		
		# Get number of iterations
		n_iter <- nrow(tmat_mat)
	
		# Get number of bodies
		n_bodies <- ncol(tmat_mat) / 16
	
		# Get body names
		body_names <- unique(gsub('_(R[0-3]{2}|[0-3]{2}|TX|TY|TZ|1)$', '', colnames(tmat_mat), ignore.case=TRUE))
		
		# Capitalize to set correct order
		colnames(tmat_mat) <- toupper(colnames(tmat_mat))
		
		# Define column order
		col_order <- c(sapply(toupper(body_names), 'paste0', c('_R11', '_R12', '_R13', '_01', '_R21', '_R22', '_R23', '_02', '_R31', '_R32', '_R33', '_03', '_TX', '_TY', '_TZ', '_1')))
		
		# Set column order
		tmat_mat <- tmat_mat[, col_order]

		# Convert to array
		tmat <- array(t(tmat_mat), dim=c(4,4,n_bodies,n_iter), dimnames=list(NULL, NULL, body_names, NULL))

		# Convert NaNs to NA
		tmat[is.na(tmat)] <- NA

	}else{
		tmat <- NULL
	}

	# Fill coordinate array
	if(has_xyz){

		# Get columns
		xyz_mat <- matrix(as.numeric(read_matrix[, xyz_cols]), nrow=nrow(read_matrix), ncol=sum(xyz_cols), dimnames=list(row_names, colnames(read_matrix)[xyz_cols]))

		# Convert XYZ matrix to array
		xyz <- mat2arr(xyz_mat)
		
		# Get number of iterations
		n_iter <- dim(xyz)[3]
	
		# Convert NaNs to NA
		xyz[is.na(xyz)] <- NA

	}else{
		xyz <- NULL
	}

	# Set
	rlist <- list()
	rlist[['tmat']] <- tmat
	rlist[['xyz']] <- xyz
	
	# Set classes
	class(rlist) <- 'motion'
	if(!is.null(tmat)) class(rlist[['tmat']]) <- 'tmat'
	if(!is.null(xyz)) class(rlist[['xyz']]) <- 'xyz'

	# Add extra info columns
	if(has_info){
		for(info_col in which(info_cols)){
		
			# Get values
			info_col_vals <- read_matrix[, info_col]
			
			#
			info_col_vals[grepl('[ ]?NA[ ]?', info_col_vals)] <- NA
			
			# Find first non-NA value
			first_nna <- which(!is.na(info_col_vals))[1]
			
			# Check if numeric
			if(length(first_nna) > 0){

				if(is.na(suppressWarnings(as.numeric(info_col_vals[first_nna])))){
					val_type <- 'character'
				}else{
					val_type <- 'numeric'
				}
				
				# Set correct class
				if(val_type == 'numeric') info_col_vals <- as(info_col_vals, val_type)

				# Check for logical
				if(val_type == 'character' && grepl('^FALSE$|^[ ]?TRUE$', info_col_vals)) val_type <- 'logical'

				# Format character class
				if(val_type == 'character'){

					# Remove default numeric names
					info_col_vals <- setNames(info_col_vals, NULL)
					
					# Remove leading and following spaces if ' , ' separator
					if(sep == ' , ') info_col_vals <- gsub('(^[ ])|([ ]$)', '', info_col_vals)
				}

				# Format logical class
				if(val_type == 'logical'){

					# Remove spaces
					info_col_vals <- gsub('[ ]', '', info_col_vals)

					# Set class
					info_col_vals <- as(info_col_vals, val_type)
				}
			}

			# Get number of iterations
			n_iter <- length(info_col_vals)

			rlist[[colnames(read_matrix)[info_col]]] <- info_col_vals
		}
	}
	
	# Set number of iterations
	rlist$n.iter <- n_iter
	
	return(rlist)

	# Get number of frames from first file
	if(file_format == 'txt'){

	}else{

		# Read csv

		# If transformations
		if(grepl('_R11$', colnames(read_matrix)[1])){

			# Check for non-transformation columns
			tm_grepl <- grepl('_(R[1-3]{2}|[0-3]{2}|1|TX|TY|TZ)', colnames(read_matrix), ignore.case=TRUE)
		
			# Get tm matrix
			tm_matrix <- read_matrix[, tm_grepl]

			# Transformation matrix
			tmat <- matrix(suppressWarnings(as.numeric(tm_matrix)), nrow(tm_matrix), ncol(tm_matrix), 
				dimnames=dimnames(tm_matrix))
			
			# Convert transformation matrix into array
			tm_arr <- tmmat2arr(tmat)

			# Set return list
			rlist <- list('tm.arr'=tm_arr)
			
			# Get non-transformation column names
			if(sum(!tm_grepl) > 0){

				# Add column as list element
				for(non_tm_colname in colnames(read_matrix)[!tm_grepl]) rlist[[non_tm_colname]] <- gsub('^[ ]*|[ ]*$', '', read_matrix[, non_tm_colname])

				return(rlist)
			}else{

				return(rlist$tm.arr)
			}

		}else if(colnames(read_matrix)[1] == 'X'){

			read_matrix <- read_matrix[, 2:ncol(read_matrix)]
			
			if('time' %in% colnames(read_matrix)){
				tmta <- time_mat_to_arr(read_matrix)
			}else{
				return(mat2arr(read_matrix))
			}

		}else if(colnames(read_matrix)[1] == 'Frame'){

			# Remove columns with all NA values
			read_matrix <- read_matrix[, colSums(!is.na(read_matrix)) > 0]

			# Remove first column if frame number
			if(colnames(read_matrix)[1] == 'Frame') read_matrix <- read_matrix[, 2:ncol(read_matrix)]

			# Convert XYZ matrix to array
			arr <- mat2arr(read_matrix, pattern='(_x|_y|_z)$')

			# Remove filter frequency from rownames, if present
			dimnames(arr)[[1]] <- gsub('_[0-9]+Hz', '', dimnames(arr)[[1]])
	
			# Sort markers alphabetically by name
			arr <- arr[sort(dimnames(arr)[[1]]), , ]

			return(arr)
		}
	}
	
	lm_array_ex <- tmta$arr

	# Set number of iterations
	num_iter <- dim(lm_array_ex)[3]
	
	# Set landmark names
	if(is.null(landmark.names)) landmark.names <- dimnames(lm_array_ex)[[1]]

	# Create array for coordinates from all strikes
	if(multiple.as == 'list' && length(file) > 1){
		marker_array <- list()
		for(ii in 1:length(file)){
			marker_array[[ii]] <- array(NA, dim=c(dim(lm_array_ex)[1:2], num_iter), dimnames=list(landmark.names, letters[24:26], NULL))
		}
	}else if(multiple.as == 'array' && length(file) > 1){

		# Get file names from paths
		file_names <- rep(NA, length(file))
		for(i in 1:length(file)){
			str_split <- strsplit(file[i], split='/')[[1]]
			file_names[i] <- gsub('[.](txt|csv)$', '', str_split[length(str_split)], ignore.case=TRUE)
		}

		# Create array
		marker_array <- array(NA, dim=c(dim(lm_array_ex)[1:2], num_iter, length(file)), dimnames=list(landmark.names, letters[24:26], NULL, file_names))

	}else{
		marker_array <- array(NA, dim=c(dim(lm_array_ex)[1:2], num_iter*length(file)), dimnames=list(landmark.names, letters[24:26], NULL))
	}

	# Read coordinates into list/array
	for(i in 1:length(file)){

		# Read coordinates
		if(file_format == 'txt'){
			tmta <- time_mat_to_arr(as.matrix(read.table(file[i])))
		}else{
			read_matrix <- as.matrix(read.csv(file[i]))
			if(colnames(read_matrix)[1] == 'X') read_matrix <- read_matrix[, 2:ncol(read_matrix)]
			tmta <- time_mat_to_arr(read_matrix)
		}

		lm_array <- tmta$arr

		# Add to array/list
		if(multiple.as == 'list' && length(file) > 1){
			marker_array[[i]][landmark.names, , ] <- lm_array[landmark.names, , ]
		}else if(multiple.as == 'array' && length(file) > 1){
			marker_array[landmark.names, , , i] <- lm_array[landmark.names, , ]
		}else{
			# Set frame indices
			frames_idx <- ((i-1)*num_iter+1):(i*num_iter)

			# Add to array
			marker_array[landmark.names, , frames_idx] <- lm_array[landmark.names, , ]
		}
	}
	
	list(
		'xyz'=marker_array, 
		'time'=tmta$times
	)
}