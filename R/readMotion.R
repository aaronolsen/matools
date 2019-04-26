readMotion <- function(file, nrows = -1, vectors.as = c('list', 'data.frame'), vectors.name = 'vectors',
	xyz.pattern = '[_|.](|X|Y|Z)$', tm.pattern = '_(R[0-3]{2}|[0-3]{2}|TX|TY|TZ|1)$'){

	## Reads in matrix of coordinates over time, with or without time column, or 
	## transformation matrices. File type is detected based on whether first column name 
	## ends in R11

	if(length(file) > 1){

		# Start with NULL motion object
		all_motion <- NULL

		# Multiple files
		for(i in 1:length(file)){

			if(!file.exists(file[i])) stop(paste0("'", file[i], "' not found."))

			# Read each file
			read_motion <- readMotion(file[i], nrows=nrows, vectors.as=vectors.as, vectors.name=vectors.name,
				xyz.pattern=xyz.pattern, tm.pattern=tm.pattern)

			# Add rows to all motion object
			all_motion <- addMotion(read_motion, all_motion)
		}
	
		return(all_motion)

	}else{
		
		if(!file.exists(file)) stop(paste0("'", file, "' not found."))

		if(file.info(file)$isdir){

			# Remove end backslash for consistency (not sure if this is necessary)
			file <- gsub('/$', '', file)

			# Get files in directory
			list_files <- list.files(file)
			
			return(readMotion(paste0(file, '/', list_files), nrows=nrows, vectors.as=vectors.as, vectors.name=vectors.name,
				xyz.pattern=xyz.pattern, tm.pattern=tm.pattern))
		}
	}

	# Set disallowed input names
	disallowed_names <- c('replace.rows', 'remove.rows', 'n.iter')

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

	# Get first line
	first_line <- readLines(file[1], n=1)
	
	# Get column names (using read.csv converts '-' into '.')
	col_names <- strsplit(first_line, split=sep)[[1]]
	if(quote != '') col_names <- gsub(quote, '', col_names)

	# Set column names
	if(length(col_names) == ncol(read_matrix) + 1) col_names <- col_names[2:length(col_names)]
	colnames(read_matrix) <- col_names

	# If first column is X, remove
	row_names <- NULL
	if(colnames(read_matrix)[1] == 'X'){
		row_names <- read_matrix[, 1]
		read_matrix <- read_matrix[, 2:ncol(read_matrix)]
	}

	# Sort column names and remove columns with no name
	col_names <- sort(colnames(read_matrix))
	col_names <- col_names[!col_names == '']

	## Correct for bug in XMALab export where there are two 'TY' columns
	# Check if there are transformations
	tmat_cols <- grepl(paste0(substr(tm.pattern, 1, nchar(tm.pattern)-1), '_[0-9]*Hz$'), colnames(read_matrix), ignore.case=TRUE)
	if(sum(tmat_cols) > 0){
		
		# Find all TY columns
		ty_cols_match <- which(grepl('_TY', colnames(read_matrix)[tmat_cols]))
		
		if(length(ty_cols_match) > 1){

			# Look for consecutive 'TY' columns
			for(i in 2:length(ty_cols_match)){
				
				# Replace second with 'TZ'
				if(ty_cols_match[i-1] == ty_cols_match[i]-1){
					colnames(read_matrix)[ty_cols_match[i]] <- gsub('_TY$', '_TZ', colnames(read_matrix)[ty_cols_match[i]])
					colnames(read_matrix)[ty_cols_match[i]] <- gsub('_TY_', '_TZ_', colnames(read_matrix)[ty_cols_match[i]])
				}
			}
		}
		
		col_names <- colnames(read_matrix)
	}

	# Sort column names
	read_matrix <- read_matrix[, col_names]

	# Check if there are transformations
	tmat_cols <- grepl(tm.pattern, colnames(read_matrix), ignore.case=TRUE)

	# Try looking for transformation columns with smoothing frequency
	if(sum(tmat_cols) == 0){
		
		# Try
		tmat_cols <- grepl(paste0(substr(tm.pattern, 1, nchar(tm.pattern)-1), '_[0-9]*Hz$'), colnames(read_matrix), ignore.case=TRUE)

		# If found, remove smoothing frequency from column names
		if(sum(tmat_cols) > 0){
			colnames(read_matrix) <- gsub('_[0-9]*Hz$', '', colnames(read_matrix))
		}
	}
	
	# Check if there are xyz coordinates
	xyz_cols <- grepl(xyz.pattern, colnames(read_matrix), ignore.case=TRUE)

	# Check if any Z columns
	z_cols <- grepl('[_|.](Z)$', colnames(read_matrix), ignore.case=TRUE)

	# Check for additional info columns
	info_cols <- tmat_cols+xyz_cols == 0

	# Set logicals
	has_tmat <- ifelse(sum(tmat_cols) > 0, TRUE, FALSE)
	has_xyz <- ifelse(sum(xyz_cols) > 0, TRUE, FALSE)
	if(has_xyz && sum(z_cols) == 0){
		has_xyz <- FALSE
		has_xy <- TRUE
	}else{
		has_xy <- FALSE
	}
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
		body_names <- unique(gsub(tm.pattern, '', colnames(tmat_mat), ignore.case=TRUE))
		
		# Capitalize to set correct order
		colnames(tmat_mat) <- toupper(colnames(tmat_mat))
		
		# Define column order
		col_order <- c(sapply(toupper(body_names), 'paste0', c('_R11', '_R12', '_R13', '_01', '_R21', '_R22', '_R23', '_02', '_R31', '_R32', '_R33', '_03', '_TX', '_TY', '_TZ', '_1')))
#		print(colnames(tmat_mat))
		#print(col_order[!col_order %in% colnames(tmat_mat)])
		
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
		xyz <- mat2arr(xyz_mat, pattern=xyz.pattern, ignore.case=TRUE)
		
		# Get number of iterations
		n_iter <- dim(xyz)[3]
	
		# Convert NaNs to NA
		xyz[is.na(xyz)] <- NA

	}else{
		xyz <- NULL
	}

	# Fill 2D coordinate array
	if(has_xy){

		# Get xy columns
		xy_cols <- grepl('[_|.](X|Y)$', colnames(read_matrix), ignore.case=TRUE)

		# Get columns
		xy_mat <- matrix(as.numeric(read_matrix[, xy_cols]), nrow=nrow(read_matrix), ncol=sum(xy_cols), dimnames=list(row_names, colnames(read_matrix)[xy_cols]))

		# Convert XYZ matrix to array
		xy <- mat2arr(xy_mat, pattern='[.|_](x|y)$', ignore.case=TRUE)

		# If 'cam#'	present, split into cams
		if(grepl('cam[0-9]+', dimnames(xy)[[1]][1])){

			# Get dimnames 1			
			dimnames1 <- dimnames(xy)[[1]]

			# Get camera names
			cam_names <- rep(NA, dim(xy)[1])
			for(i in 1:dim(xy)[1]){
				reg_expr <- regexpr('cam[0-9]+', dimnames1[i])
				cam_names[i] <- substr(dimnames1[i], reg_expr, reg_expr+attr(reg_expr, 'match.length'))
			}
			cam_names_unique <- sort(unique(cam_names))

			# Create point names without cameras
			dimnames1_wocam <- gsub('_cam[0-9]+$', '', dimnames1)
			dimnames1_wocam <- gsub('_cam[0-9]+_', '_', dimnames1_wocam)
			dimnames1_wocam_unique <- unique(dimnames1_wocam)
			
			# Create new array with cameras
			xy_cam <- array(NA, dim=c(length(dimnames1_wocam_unique), dim(xy)[2],length(cam_names_unique),dim(xy)[3]), 
				dimnames=list(dimnames1_wocam_unique, NULL, cam_names_unique, NULL))
			
			# Fill array
			for(cam_name in cam_names_unique) xy_cam[dimnames1_wocam[cam_names == cam_name], , cam_name, ] <- xy[cam_names == cam_name, , ]

			# Replace array			
			xy <- xy_cam

			# Get number of iterations
			n_iter <- dim(xy)[4]
		}else{

			# Get number of iterations
			n_iter <- dim(xy)[3]
		}

		# Convert NaNs to NA
		xy[is.na(xy)] <- NA

	}else{
		xy <- NULL
	}

	# Set
	rlist <- list()
	rlist[['tmat']] <- tmat
	rlist[['xyz']] <- xyz
	rlist[['xy']] <- xy
	
	# Set classes
	class(rlist) <- 'motion'
	if(!is.null(tmat)) class(rlist[['tmat']]) <- 'tmat'
	if(!is.null(xyz)) class(rlist[['xyz']]) <- 'xyz'
	if(!is.null(xy)) class(rlist[['xy']]) <- 'xy'

	# Add extra info columns
	if(has_info){
	
		val_types <- rep(NA, sum(info_cols))
		info_cols_idx <- which(info_cols)

		for(i in 1:sum(info_cols)){
		
			#
			info_col <- info_cols_idx[i]

			#
			if(colnames(read_matrix)[info_col] %in% disallowed_names) stop(paste0("Columns in the motion file cannot have any of the following names since they are used for internal operations: ", paste0(disallowed_names, collapse=', ')))
		
			# Get values
			info_col_vals <- read_matrix[, info_col]
			
			#
			info_col_vals[grepl('[ ]?NA[ ]?', info_col_vals)] <- NA
			
			# Find first non-NA value
			first_nna <- which(!is.na(info_col_vals))[1]
			
			# Check if numeric
			if(length(first_nna) > 0){

				if(is.na(suppressWarnings(as.numeric(info_col_vals[first_nna])))){
					val_types[i] <- 'character'
				}else{
					val_types[i] <- 'numeric'
				}
				
				# Set correct class
				if(val_types[i] == 'numeric') info_col_vals <- as(info_col_vals, val_types[i])

				# Check for logical
				if(val_types[i] == 'character' && grepl('^FALSE$|^[ ]?TRUE$', info_col_vals)) val_types[i] <- 'logical'

				# Format character class
				if(val_types[i] == 'character'){

					# Remove default numeric names
					info_col_vals <- setNames(info_col_vals, NULL)
					
					# Remove leading and following spaces if ' , ' separator
					if(sep == ' , ') info_col_vals <- gsub('(^[ ])|([ ]$)', '', info_col_vals)
				}

				# Format logical class
				if(val_types[i] == 'logical'){

					# Remove spaces
					info_col_vals <- gsub('[ ]', '', info_col_vals)

					# Set class
					info_col_vals <- as(info_col_vals, val_types[i])
				}
			}

			# Get number of iterations
			n_iter <- length(info_col_vals)

			rlist[[colnames(read_matrix)[info_col]]] <- info_col_vals
		}
		
		if(vectors.as[1] == 'data.frame'){

			rlist_to_df <- do.call(cbind.data.frame, rlist[info_cols])
			rlist[info_cols] <- NULL
			rlist[[vectors.name]] <- rlist_to_df

			for(i in 1:ncol(rlist[[vectors.name]])){
			
				# Get column name
				col_name <- colnames(rlist[[vectors.name]])[i]
	
				# If character, convert from factor to character
				if(val_types[i] == 'character') rlist[[vectors.name]][[col_name]] <- as.character(rlist[[vectors.name]][[col_name]])
			}
		}
	}
	
	# Set number of iterations
	rlist$n.iter <- n_iter
	
	return(rlist)
}

print.motion <- function(x){
	
	rc <- ''
	
	# Get names of objects
	x_names <- names(x)

	# Print transformation matrix details	
	if('n.iter' %in% x_names){
		rc <- c(rc, paste0('Number of iterations: ', x$n.iter, '\n'))
	}
	if('tmat' %in% x_names){
		rc <- c(rc, paste0('$tmat (', paste0(dim(x$tmat), collapse='x'), ')', '\n'))
		if(length(dim(x$tmat)) > 2 && !is.null(dimnames(x$tmat)[[3]])){
			dots <- ''
			if(dim(x$tmat)[3] > 10) dots <- paste0('\n\t... and ', dim(x$tmat)[3]-10, ' more bodies')
			rc <- c(rc, paste0('\t', paste0(dimnames(x$tmat)[[3]][1:min(10,dim(x$tmat)[3])], collapse='\n\t'), dots, '\n'))
		}
	}

	# Print transformation matrix details	
	if('xyz' %in% x_names){
		rc <- c(rc, paste0('$xyz (', paste0(dim(x$xyz), collapse='x'), ')', '\n'))
		if(!is.null(dimnames(x$xyz)[[1]])){
			dots <- ''
			if(dim(x$xyz)[1] > 10) dots <- paste0('\n\t... and ', dim(x$xyz)[1]-10, ' more rows')
			rc <- c(rc, paste0('\t', paste0(dimnames(x$xyz)[[1]][1:min(10,dim(x$xyz)[1])], collapse='\n\t'), dots, '\n'))
		}
	}

	# Print transformation matrix details	
	if('xy' %in% x_names){
		rc <- c(rc, paste0('$xy (', paste0(dim(x$xy), collapse='x'), ')', '\n'))
		if(!is.null(dimnames(x$xy)[[1]])) rc <- c(rc, paste0('\t', paste0(dimnames(x$xy)[[1]], collapse='\n\t'), '\n'))
	}

	info_names <- x_names[!x_names %in% c('xyz', 'xy', 'tmat', 'n.iter', 'replace.rows', 'remove.rows')]

	if(length(info_names) > 0){
		xlist_to_df <- do.call(cbind.data.frame, x[info_names])
		colnames(xlist_to_df) <- paste0('$', colnames(xlist_to_df))
		#rc <- c(rc, 'Other objects\n')
		rc <- c(rc, paste0(paste0(capture.output(print(head(xlist_to_df))), collapse='\n'), '\n'))
	}

	cat(rc, sep='')
}