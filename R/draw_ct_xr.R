draw_ct_xr <- function(file, ct_arr, xr_arr, cs_arr, path.connect, duration = 10){

	# Get xyz limits for CT and X-ray markers
	ranges <- rbind(apply(ct_arr, 2, 'range', na.rm=TRUE), apply(xr_arr, 2, 'range', na.rm=TRUE))

	# Remove frames where all CT markers are NA
	xr_arr <- xr_arr[, , colSums(is.na(ct_arr[,1,])) < dim(ct_arr)[1]]
	ct_arr <- ct_arr[, , colSums(is.na(ct_arr[,1,])) < dim(ct_arr)[1]]

	# Read path connect into list
	read_lines <- readLines(path.connect)
	plist_names <- list()
	for(line_num in 1:length(read_lines)) plist_names[[length(plist_names)+1]] <- strsplit(read_lines[line_num], ',[ ]?')[[1]]

	# Save CT marker names for reference
	ct_marker_names <- dimnames(ct_arr)[[1]]
	xr_marker_names <- dimnames(xr_arr)[[1]]

	# Create path connect list
	plist_ct <- list()
	plist_xr <- list()
	for(i in 1:length(plist_names)){
		v_ct <- c()
		v_xr <- c()
		for(j in c(1:length(plist_names[[i]]))){
			v_ct <- c(v_ct, which(ct_marker_names == plist_names[[i]][j]))
			v_xr <- c(v_xr, which(xr_marker_names == plist_names[[i]][j]))
		}
		plist_ct[[i]] <- v_ct
		plist_xr[[i]] <- v_xr
	}

	# Identify virtual markers
	is_virtual <- setNames(rep(FALSE, dim(xr_arr)[1]), dimnames(xr_arr)[[1]])
	is_virtual[grepl('-', dimnames(xr_arr)[[1]])] <- TRUE

	# Create new svg file
	svg.new(file=file, animate.duration=duration)

	# Draw CT markers
	svg.points(ct_arr)

	# Draw coordinate systems
	if(!is.null(cs_arr)){

		# Get XYZ centroid size
		xr_xyz <- mean(dppt(xr_arr[, , 1], colMeans(xr_arr[, , 1], na.rm=TRUE)), na.rm=TRUE)

		# Set coordinate system arrow properties
		arrow_len <- xr_xyz*0.25
		arrowhead_len <- xr_xyz*0.03

		# Set coordinate system colors
		col_cs <- c('red', 'green', 'blue')

		for(cs_name in dimnames(cs_arr)[[3]]){

			# Each arrow
			for(j in 2:4){

				# Scale vectors based on plot size so they are visible in animation
				cs_arr_sub <- cs_arr[c(1,j), , cs_name, ]
				cs_arr_sub[2, , ] <- cs_arr_sub[1, , ] + (cs_arr_sub[2, , ]-cs_arr_sub[1, , ])*arrow_len
			
				# Draw coordinate systems as arrows
				#if(cs_name == 'Urohyal' && j == 2) print(cs_arr_sub)
				svg.arrows(x=cs_arr_sub, len=arrowhead_len, 
					col=col_cs[j-1], layer='Coordinate systems', z.index=1)
			}
		}
	}

	# Connect CT markers
	svg.pathsC(plist_ct)
	svg.pathsC(plist_xr, col.stroke='red')

	# Draw X-ray markers
	if(sum(is_virtual) > 0)	svg.points(array(xr_arr[is_virtual, , ], dim=c(sum(is_virtual), dim(xr_arr)[2:3])), cex=1, col='orange')
	if(sum(!is_virtual) > 0) svg.points(xr_arr[!is_virtual, , ], cex=1, col='red')

	# Create frame
	svg.frame(ranges, z.index=-1)

	svg.close()
	
}