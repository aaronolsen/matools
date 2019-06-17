unifyMotion <- function(motion, xyz.mat, unify.spec, regexp = FALSE, 
	print.progress = TRUE, print.progress.iter = c(1), verbose = FALSE,
	replace.xyz = TRUE, plot.diag = NULL, cp.use = TRUE, cp.axis.with.vp = FALSE){

	# Set point array
	if(is.list(motion)){
		input_list <- TRUE
		xr_arr <- motion$xyz
		n_iter <- motion$n.iter
		frames <- 1:n_iter
		if(!is.null(motion$frame)) frames <- motion$frame
		if(!is.null(motion$Frame)) frames <- motion$Frame
	}else{
		input_list <- FALSE
		xr_arr <- xyz
		n_iter <- dim(xyz)[3]
		frames <- 1:n_iter
	}
	
	# Set matrix coordinates
	ct_mat <- xyz.mat

	# Remove CT coordinates that are NA
	ct_mat <- ct_mat[!is.na(ct_mat[, 1]), ]
	
	# Set CT markers
	ct_markers <- rownames(ct_mat)

	# Parse unify specifications and return in unification order
	ulist <- parse_unify_spec(ulist=unify.spec, xr_arr=xr_arr, ct_mat=ct_mat, regexp=regexp, print.progress=print.progress)

	# Get element names from list
	body_names <- names(ulist)

	# Print ulist
	if(print.progress && verbose){

		cat(paste0('Parse unify specifications\n'))

		for(elem_name in names(ulist)){
			cat(paste0('\t', elem_name, '\n'))
			for(l_type in c('Align', 'Point', 'Plane', 'Transform')){
				if(is.null(ulist[[elem_name]][[tolower(l_type)]])) next
				if(l_type == 'Plane'){
					other_name <- names(ulist[[elem_name]][['plane']])[names(ulist[[elem_name]][['plane']]) != elem_name]
					cat(paste0('\t\t', l_type, '\n'))
					cat(paste0('\t\t\tMarkers to move into plane (', elem_name, '): ', paste0(unlist(ulist[[elem_name]][[tolower(l_type)]][[elem_name]]), collapse=', '), '\n'))
					cat(paste0('\t\t\tMarkers used to define plane (', other_name, '): ', paste0(unlist(ulist[[elem_name]][[tolower(l_type)]][[other_name]]), collapse=', '), '\n'))
				}else{
					cat(paste0('\t\t', l_type, ': ', paste0(unlist(ulist[[elem_name]][[tolower(l_type)]]), collapse=', '), '\n'))
				}
			}
		}
	}

	# Create transformation matrix array
	tm_arr <- array(NA, dim=c(4, 4, length(body_names), dim(xr_arr)[3]), dimnames=list(NULL, NULL, body_names, NULL))

	# Points in motion markers
	ct_arr_pts <- ct_markers[ct_markers %in% dimnames(xr_arr)[[1]]]
	
	# And point-to markers
	point_to_pts <- c()
	for(body_name in names(ulist)){
		if(!is.null(ulist[[body_name]][['transform']])) point_to_pts <- c(point_to_pts, ulist[[body_name]][['transform']])
	}
	ct_arr_pts <- c(ct_arr_pts, point_to_pts)
	ct_arr_pts <- sort(unique(ct_arr_pts))

	# If non NULL create array to store results of constraint plane transformations
	unify_cp <- NULL
	if(cp.use && !is.null(unify_cp)){
		col_names <- c('pre.min', 'pre.mean', 'pre.max', 'angle', 'post.min', 'post.mean', 'post.max')
		cp_array <- array(NA, dim=c(length(unify_cp$id), length(col_names), n_iter), dimnames=list(NULL, col_names, NULL))
	}

	# Add virtual markers to xr_arr
	xr_arr_n <- array(NA, dim=c(dim(xr_arr)[1]+length(ct_arr_pts), dim(xr_arr)[2], dim(xr_arr)[3]), 
		dimnames=list(c(dimnames(xr_arr)[[1]], ct_arr_pts), dimnames(xr_arr)[[2]], dimnames(xr_arr)[[3]]))
	xr_arr_n[dimnames(xr_arr)[[1]], , ] <- xr_arr

	# Create unification error matrix
	errors <- matrix(NA, nrow=dim(xr_arr_n)[3], ncol=length(body_names), dimnames=list(NULL, body_names))

	## Iterate through each frame
	min_iter <- 1
	max_iter <- dim(xr_arr_n)[3]
	for(iter in min_iter:max_iter){
	
		# Set print progress for iteration
		if(print.progress && iter %in% print.progress.iter){ print_progress <- TRUE }else{ print_progress <- FALSE }

		# Print progress
		if(print_progress) cat(paste0('Iteration: ', iter, '\n'))
		
		# Unify markers
		for(body_name in body_names){

			# Get non-NA motion markers
			mo_markers <- dimnames(xr_arr_n)[[1]][!is.na(xr_arr_n[, 1, iter])]

			#if(!body_name %in% c('SuspensoriumL', 'Neurocranium')) next

			if(print_progress) cat(paste0('\t', body_name, '\n'))

			# Find align markers in motion array
			align_markers <- ulist[[body_name]][['align']][ulist[[body_name]][['align']] %in% mo_markers]
			
			# Skip if no align markers
			if(length(align_markers) == 0) next

			# Get point markers
			if(!is.null(ulist[[body_name]][['point']])){

				# Find point-to markers in both motion and CT
				point_markers <- unlist(ulist[[body_name]][['point']])
				point_markers <- point_markers[point_markers %in% mo_markers]

			}else{
				point_markers <- c()
			}
			
			# Get plane markers
			if(!is.null(ulist[[body_name]][['plane']])){
			
				# Get name of plane body
				plane_body <- names(ulist[[body_name]][['plane']])[names(ulist[[body_name]][['plane']]) != body_name]

				# Get plane markers
				plane_markers <- unlist(ulist[[body_name]][['plane']][[plane_body]])
				plane_markers <- plane_markers[plane_markers %in% mo_markers]

				# Fit plane to points
				fit_plane <- fitPlane(xr_arr_n[plane_markers,, iter])

				# Get markers to move into plane
				markers_in_plane <- ulist[[body_name]][['plane']][[body_name]]

				# Set plane marker indices
				plane_marker_idx <- (length(point_markers)+1):(length(point_markers)+length(markers_in_plane))

			}else{
				plane_markers <- c()
				markers_in_plane <- c()
			}
			
			# Get unique list of all markers
			all_markers <- unique(c(align_markers, point_markers, markers_in_plane))

			if(length(align_markers) >= 2){

				# Print progress
				if(print_progress) cat(paste0('\t\tAlign CT markers using ', length(align_markers), ' motion markers: "', paste0(align_markers, collapse='", "'), '"\n'))

				# Transform CT markers to correspond with motion markers
				align <- bestAlign(xr_arr_n[align_markers,, iter], ct_mat[unique(c(align_markers, point_markers, markers_in_plane)),], sign=1)
			}

			if(length(point_markers) == 0 && (length(plane_markers) == 0 || length(markers_in_plane) == 0)){
			
				## Only use align markers
				# Save transformation matrix
				tm_arr[, , body_name, iter] <- align$tmat

			}else{

				# Single alignment marker
				if(length(align_markers) == 1){

					# Create transformed CT mat
					ct_mat_t <- ct_mat[c(align_markers, point_markers, markers_in_plane), , drop=FALSE]

					# Set real marker as center
					center <- xr_arr_n[align_markers,, iter]

					# Get translation vector based on single real marker
					translate_tmat <- diag(4)
					translate_tmat[1:3, 4] <- center - ct_mat[align_markers, ]

					# Translate body points before rotating
					ct_mat_t <- ct_mat_t + matrix(translate_tmat[1:3, 4], nrow(ct_mat_t), 3, byrow=TRUE)

					if(!is.null(ulist[[body_name]][['plane']])){
					
						# Print progress
						if(print_progress){
							cat(paste0('\t\tAlign CT markers using 1 motion marker: "', align_markers, '"\n'))
							if(length(point_markers) > 0) cat(paste0('\t\tOptimize rotation using point-to markers: "', paste0(point_markers, collapse='", "'), '"\n'))
							cat(paste0('\t\tPlacing marker(s) "', paste0(markers_in_plane, collapse='", "'), '" in plane...\n\t\t\tdefined by markers: "', paste0(plane_markers, collapse='", "'), '"\n'))
						}
						
						# Find initial error
						rotate_error_init <- ref_rotate_error(c(0,0,0), center=center,
							fit.plane.point=fit_plane$Q, fit.plane.normal=fit_plane$N, plane.marker.idx=plane_marker_idx,
							ref.points=ct_mat_t[c(point_markers, markers_in_plane),], fit.points=xr_arr_n[point_markers,, iter])
						
						# Optimize points by rotating 3 axes about real marker
						rotation_fit <- tryCatch(
							expr={
								nlminb(start=c(0,0,0), objective=ref_rotate_error, lower=-2*pi, upper=2*pi, center=center, 
									fit.plane.point=fit_plane$Q, fit.plane.normal=fit_plane$N, plane.marker.idx=plane_marker_idx, 
									ref.points=ct_mat_t[c(point_markers, markers_in_plane),], fit.points=xr_arr_n[point_markers,, iter])
							},
							error=function(cond) {print(cond);return(NULL)},
							warning=function(cond) {print(cond);return(NULL)}
						)

					}else{

						# Print progress
						if(print_progress){
							cat(paste0('\t\tAlign CT markers using 1 motion marker: "', align_markers, '"\n'))
							cat(paste0('\t\tOptimize rotation using point markers: "', paste0(point_markers, collapse='", "'), '"\n'))
						}

						# Find initial error
						rotate_error_init <- ref_rotate_error(c(0,0,0), center=center,  
							ref.points=ct_mat_t[point_markers,], fit.points=xr_arr_n[point_markers,, iter])

						# Optimize points by rotating 3 axes about real marker
						rotation_fit <- tryCatch(
							expr={
								nlminb(start=c(0,0,0), objective=ref_rotate_error, lower=-2*pi, upper=2*pi, center=center, 
									ref.points=ct_mat_t[point_markers,], fit.points=xr_arr_n[point_markers,, iter])
							},
							error=function(cond) {print(cond);return(NULL)},
							warning=function(cond) {print(cond);return(NULL)}
						)
					}

					if(print_progress) cat(paste0('\t\tError prior to optimization: ', rotate_error_init, '; Error after optimization: ', rotation_fit$objective, '\n'))

					# Save transformation matrix using optimized angle, including initial transformation
					tmat1 <- tmat2 <- tmat3 <- diag(4)
					tmat1[1:3,4] <- center
					tmat2[1:3,1:3] <- rotationMatrixZYX_ma(rotation_fit$par)
					tmat3[1:3,4] <- -center
					tm_arr[, , body_name, iter] <- tmat1 %*% tmat2 %*% tmat3 %*% translate_tmat

				}else{

					# Print progress
					if(print_progress){
						cat(paste0('\t\tOptimize rotation about the axis defined by align markers: "', paste0(align_markers, collapse='", "'), '"\n'))
						if(length(point_markers) > 0) cat(paste0('\t\tOptimize rotation using point-to marker(s): "', paste0(point_markers, collapse='", "'), '"\n'))
					}

					# Get axis for refining rotation
					if(length(align_markers) == 2){

						# Set axis and center of rotation
						raxis <- uvector_ma(xr_arr_n[align_markers[2],, iter] - xr_arr_n[align_markers[1],, iter])
						center <- xr_arr_n[align_markers[2],, iter]

					}else{

						# Print progress
						#if(print_progress) cat(paste0('\t\tOptimize rotation about axis fit to align markers to point to marker(s): "', paste0(point_markers, collapse='", "'), '"\n'))

						# Set axis and center of rotation
						fit_line <- fitLine3D_ma(align$mat[align_markers, ])
						raxis <- uvector_ma(fit_line$p2-fit_line$p1)
						center <- fit_line$p1
					}

					# Create transformed CT mat
					ct_mat_t <- align$mat

					# Find initial error
					if(!is.null(ulist[[body_name]][['plane']])){

						# Print progress
						if(print_progress) cat(paste0('\t\tPlacing marker(s) "', paste0(markers_in_plane, collapse='", "'), '" in plane...\n\t\t\tdefined by markers: "', paste0(plane_markers, collapse='", "'), '"\n'))

						rotate_error_init <- ref_rotate_error(0, center=center, axis=raxis, 
							fit.plane.point=fit_plane$Q, fit.plane.normal=fit_plane$N, plane.marker.idx=plane_marker_idx, 
							ref.points=ct_mat_t[c(point_markers, markers_in_plane),], fit.points=xr_arr_n[point_markers,, iter])
						
						# Run optimization
						rotation_fit <- tryCatch(
							expr={
								nlminb(start=0, objective=ref_rotate_error, lower=-2*pi, upper=2*pi, center=center, 
									fit.plane.point=fit_plane$Q, fit.plane.normal=fit_plane$N, plane.marker.idx=plane_marker_idx, 
									axis=raxis, ref.points=ct_mat_t[c(point_markers, markers_in_plane),], fit.points=xr_arr_n[point_markers,, iter])
							},
							error=function(cond) {print(cond);return(NULL)},
							warning=function(cond) {print(cond);return(NULL)}
						)

					}else{

						rotate_error_init <- ref_rotate_error(0, center=center, axis=raxis, 
							ref.points=ct_mat_t[point_markers,], fit.points=xr_arr_n[point_markers,, iter])

						# Run optimization
						rotation_fit <- tryCatch(
							expr={
								nlminb(start=0, objective=ref_rotate_error, lower=-2*pi, upper=2*pi, center=center, 
									axis=raxis, ref.points=ct_mat_t[point_markers,], fit.points=xr_arr_n[point_markers,, iter])
							},
							error=function(cond) {print(cond);return(NULL)},
							warning=function(cond) {print(cond);return(NULL)}
						)
					}

					if(print_progress) cat(paste0('\t\tError prior to optimization: ', rotate_error_init, '; Error after optimization: ', rotation_fit$objective, '\n'))

					# Save transformation matrix using optimized angle, including initial transformation
					tmat1 <- tmat2 <- tmat3 <- diag(4)
					tmat1[1:3,4] <- center
					tmat2[1:3,1:3] <- tMatrixEP_ma(raxis, rotation_fit$par)
					tmat3[1:3,4] <- -center
					tm_arr[, , body_name, iter] <- tmat1 %*% tmat2 %*% tmat3 %*% align$tmat
				}
			}

			# Transform associated points
			if(!is.null(ulist[[body_name]][['transform']])){
			
				# Transform markers
				transform_markers <- ulist[[body_name]][['transform']]
				
				if(print_progress){
					cat(paste0('\t\tTransform markers: "', paste0(transform_markers, collapse='", "'), '"\n'))
				}

				# Apply transformation
				if(!is.na(tm_arr[1, 1, body_name, iter])){
					xr_arr_n[transform_markers, , iter] <- applyTransform(to=ct_mat[transform_markers, ], tmat=tm_arr[, , body_name, iter])
				}
			}


			# Find error - after optimizations
			align_markers_t <- applyTransform(to=ct_mat[align_markers, ], tmat=tm_arr[, , body_name, iter])
			if(length(align_markers) == 1){
				errors_by_row <- sqrt(sum((xr_arr_n[align_markers,, iter] - align_markers_t)^2))
				errors[iter, body_name] <- errors_by_row
			}else{
				errors_by_row <- sqrt(rowSums((xr_arr_n[align_markers,, iter] - align_markers_t)^2))
				errors[iter, body_name] <- mean(errors_by_row)
			}
			
			if(print_progress){
				cat(paste0('\t\tErrors:\n\t\t\t', paste0(align_markers, ': ', c(round(errors_by_row, 3)), collapse='\n\t\t\t'), '\n'))
			}

			## Impose any constraint planes
			# Check if body is child in a constraint plane set
			if(cp.use && body_name %in% unify_cp$child){

				# Get CP IDs
				cp_ids <- unify_cp$id[body_name == unify_cp$child]
				
				# Get matching index
				cp_match <- unify_cp$id == cp_ids[1]

				if(print_progress){
					if(length(cp_ids) == 1){
						cat(paste0('\t\tRefining motion using constraint plane "', cp_ids,'" and parent body "', unify_cp$parent[cp_match], '"\n'))
					}else{
						cat(paste0('\t\tRefining motion using constraint planes "', paste0(cp_ids, collapse=','), '", parent body "', unify_cp$parent[cp_match], '"\n'))
					}
				}
				
				# Find transformation from constraint plane and points
				cpt <- suppressWarnings(constraintPlaneTransform(ct_arr[, , iter], cp_ids, unify_cp$child.type[cp_match], axis.with.vp=cp.axis.with.vp))
				tmat <- cpt$tmat
				cp_array[cp_ids,,iter] <- cpt$dist
				#print(iter)
				#print(cpt$dist)

				if(print_progress){
					cat(paste0('\t\tApplying transformation to:\n\t\t\t', paste0(rownames(ct_mat_sub_t), collapse='\n\t\t\t'), '\n'))
				}

				# Apply transformation to points and transformation
				ct_arr[rownames(ct_mat_sub_t), , iter] <- applyTransform(ct_arr[rownames(ct_mat_sub_t), , iter], tmat)
				ct_mat_sub_t <- applyTransform(ct_mat_sub_t, tmat)
				tm_arr[, , body_name, iter] <- tmat %*% tm_arr[, , body_name, iter]

				# Get new errors
				if(unify_mode[body_name] == 1){
					beads_and_vp <- rownames(ct_mat_sub_t)[grepl('_bead', rownames(ct_mat_sub_t))]
				}else{
					beads_and_vp <- rownames(ct_mat_sub_t)[grepl('_bead|-', rownames(ct_mat_sub_t))]
				}

				dist_errors <- dppt(xr_arr_n[beads_and_vp, , iter], ct_mat_sub_t[beads_and_vp, ])
				if(print_progress){
					cat(paste0('\t\tErrors:\n\t\t\t', paste0(names(dist_errors), ': ', c(round(dist_errors, 3)), collapse='\n\t\t\t'), '\n'))
				}

				# Save new errors
				errors[iter, body_name] <- mean(dist_errors, na.rm=TRUE)
			}
		}
	}

	if(cp.use && !is.null(unify_cp) && print.progress){
		cat(paste0('\nResults of constraint plane transformations:\n'))

		for(i in 1:dim(cp_array)[1]){
		
			cat(paste0('\t', i, ') ', unify_cp$parent[i]), '(')

			if(unify_cp$parent.type[i] == 'v'){ cat('plane') }else{ cat('points') }
			cat(paste0(') - ', unify_cp$child[i], ' ('))
			if(unify_cp$child.type[i] == 'v'){ cat('plane') }else{ cat('points') }
			cat(paste0(')\n'))

			# Skip if all values are NA
			if(sum(!is.na(cp_array[i,,])) == 0){
				cat(paste0('\t\tAll values NA\n'))
				next
			}

			cat(paste0('\t\tPoint-to-plane distance from: Min:', round(min(cp_array[i,'pre.min',], na.rm=TRUE), 3), ', Mean:', round(mean(cp_array[i,'pre.mean',], na.rm=TRUE), 3), ', Max:', round(max(cp_array[i,'pre.max',], na.rm=TRUE), 3), '\n'))
			if(sum(!is.na(cp_array[i,'post.min',])) > 0){
				cat(paste0('\t\tTo: Min:', round(min(cp_array[i,'post.min',], na.rm=TRUE), 3), ', Mean:', round(mean(cp_array[i,'post.mean',], na.rm=TRUE), 3), ', Max:', round(max(cp_array[i,'post.max',], na.rm=TRUE), 3), '\n'))
				cat(paste0('\t\tAngle (deg): Min:', round(min(cp_array[i,'angle',], na.rm=TRUE), 3), ', Mean:', round(mean(cp_array[i,'angle',], na.rm=TRUE), 3), ', Max:', round(max(cp_array[i,'angle',], na.rm=TRUE), 3), '\n'))
			}
			cat(paste0('\t\tNumber of frames corrected: ', sum(!is.na(cp_array[i,'post.min',])), ' of ', n_iter, '\n'))
		}
	}

	#
	if(!is.null(plot.diag)){

		# If plot filename doesn't end in pdf, add
		if(!grepl('[.]pdf$', plot.diag, ignore.case=TRUE)) plot.diag <- paste0(plot.diag, '.pdf')

		# Set which columns to plot - not columns that are all NA
		cols_plot <- which(colSums(is.na(errors)) < nrow(errors))

		## Create error diagnostic plot
		pdf(plot.diag, height=length(cols_plot)*2.5, width=max(nrow(errors)/90, 7))
		layout(cbind(1:length(cols_plot)))
		par(mar=c(4.5,4.5,3,1))
	
		for(i in cols_plot){

			if(!any(!is.na(errors[, i]))){
				plot(c(0,1), c(0,1), type='n', main=colnames(errors)[i], 
					xlab='Frame', ylab='Unification error (mm)')
				next
			}

			plot(c(frames[1], tail(frames, 1)), c(0, max(errors[, i], na.rm=TRUE)), type='n', 
				main=paste0(colnames(errors)[i], ' (Mean: ', round(mean(errors[, i], na.rm=TRUE), 3), '; SD: ', round(sd(errors[, i], na.rm=TRUE), 3), ')'), 
				xlab='Frame', ylab='Unification error (mm)')

			abline(h=0, lty=2, col=gray(0.5))

			points(frames, errors[, i], type='l')
		}
	
		dev.off()
	
	}
	
	class(errors) <- 'unify_errors'

	if(input_list){

		motion[['xyz']] <- xr_arr_n[ct_arr_pts,,]
		motion[['tmat']] <- tm_arr

		rlist <- list(
			'motion'=motion,
			'error'=errors
		)
	}else{

		rlist <- list(
			'motion'=list('xyz'=xr_arr_n[ct_arr_pts,,], 'tmat'=tm_arr),
			'error'=errors
		)
	}

	rlist
}

print.unify_errors <- function(x, n=5){
        
	rc <- ''

	class(x) <- NULL

	# Sort by column names
	x <- x[, sort(colnames(x))]

	# Name rownames so that subset of rownames can be re-named
	rownames(x) <- 1:nrow(x)

	# Set number of rows to print
	nrows_print <- min(n, nrow(x))

	# Add spaces to rownames so that columns line up between raw values and stats
	rownames(x)[1:nrows_print] <- paste0(1:nrows_print, paste0(rep(' ', max(4-nchar(n), 0)), collapse=''))

	# Convert to dataframe
	xlist_to_df <- as.data.frame(x)
	colnames(xlist_to_df) <- paste0('$', colnames(xlist_to_df))
	rc <- c(rc, paste0(paste0(capture.output(print(head(xlist_to_df, n))), collapse='\n'), '\n'))

	if(nrow(x) > n) rc <- c(rc, paste0('... and ', nrow(x)-nrows_print, ' more rows\n'))

	# Add min, max, mean
	x_stats <- matrix(NA, 3, ncol(x), dimnames=list(c('min', 'max', 'mean'), colnames(x)))
	for(i in 1:ncol(x)){
		if(!any(!is.na(x[, i]))) next
		x_stats['min', i] <- min(x[, i], na.rm=TRUE)
		x_stats['max', i] <- max(x[, i], na.rm=TRUE)
		x_stats['mean', i] <- mean(x[, i], na.rm=TRUE)
	}

	xlist_to_df <- as.data.frame(x_stats)
	colnames(xlist_to_df) <- paste0('$', colnames(xlist_to_df))
	rc <- c(rc, paste0(paste0(capture.output(print(xlist_to_df[c('min', 'max', 'mean'),])), collapse='\n'), '\n'))

	cat(rc, sep='')
}