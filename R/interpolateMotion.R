interpolateMotion <- function(motion, fit.win.max = 50, fit.win.min = 1, int.win.max = 70, span.factor = 18, plot.diag = NULL){

	## Sort to fill smallest windows first
	## Then use interpolated values in interpolation

	# Get xyz coordinates
	input_coordinates <- FALSE
	if(is.list(motion)){	
		if(!is.null(motion$xyz)){ xyz <- motion$xyz }else{ xyz <- NULL }
		
	}else{
		input_coordinates <- TRUE
		xyz <- motion
	}

	## Interpolate coordinates
	if(!is.null(xyz)){

		# Create arrays
		xyz_inter <- xyz
		
		# Time vector
		times <- 1:dim(xyz)[3]
		
		# Vector to save whether points were interpolated
		point_has_inter <- rep(FALSE, dim(xyz)[1])
		
		# Create list to save interpolated indices
		inter_indices <- list()

		# For each point
		for(i in 1:dim(xyz)[1]){

			# Skip if all values are NA
			if(!any(!is.na(xyz[i, 1, ]))) next

			# Skip if all values are not NA
			if(!any(is.na(xyz[i, 1, ]))) next
			
			# Find NA "windows"
			na_wins <- find_NA_windows(xyz[i, 1, ], sort=TRUE)
			
			# If no windows, skip
			if(is.null(na_wins)) next
			
			# If there are at least 2 windows with one non-NA value between consecutive windows
			fit_single_poly <- setNames(rep(FALSE, dim(xyz)[1]), dimnames(xyz)[[1]])
			if(nrow(na_wins) > 2){
				
				# Get length of non-NA values between windows
				non_NA_val_len <- na_wins[2:nrow(na_wins),1] - na_wins[1:(nrow(na_wins)-1),2] - 1
				
				# If this is consistent across all windows
				if(non_NA_val_len[1] == 1 && sd(non_NA_val_len) == 0) fit_single_poly[i] <- TRUE
			}

			inter_indices[[i]] <- list()

			if(fit_single_poly[i]){

				# Set that point has some interpolation
				point_has_inter[i] <- TRUE

				# Get indices of non-NA values
				non_NA_idx <- which(!is.na(xyz[i, 1, ]))
				
				# Length of non NA
				non_NA_len <- tail(non_NA_idx, 1)-non_NA_idx[1]

				# Get indices of NA values
				NA_idx <- which(is.na(xyz[i, 1, ]))
				
				# Make sure indices don't exceed non-NA indices (ie no extrapolation beyond sequence)
				NA_idx <- NA_idx[NA_idx < tail(non_NA_idx, 1)]
				NA_idx <- NA_idx[NA_idx > non_NA_idx[1]]

				# Set span
				span <- round((tail(span.factor,1)/(0.5*non_NA_len)) + 0.01, 3)

				# For each dimension
				for(j in 1:dim(xyz)[2]){
				
					# Create dataframe
					data_frame <- data.frame(y=xyz[i, j, non_NA_idx[1]:tail(non_NA_idx, 1)], x=non_NA_idx[1]:tail(non_NA_idx, 1))
					
					# Fit loess (higher span value corresponds to lower parameter fit)
					loess_main <- suppressWarnings(loess(y ~ x, data=data_frame, span=span))
				
					# Create smoothed points
					smoothed <- rep(NA, non_NA_len)
					smoothed[non_NA_idx[1]:tail(non_NA_idx, 1)] <- predict(loess_main, data_frame)
					
					# Save which indices were interpolated
					if(j == 1) inter_indices[[i]][[1]] <- NA_idx
				
					# Fill in smoothed values
					xyz_inter[i, j, NA_idx] <- smoothed[NA_idx]
				}

			}else{

				#print(na_wins)

				# For each window
				for(win_num in 1:nrow(na_wins)){
			
					# Check that window is narrower than max
					if((na_wins[win_num, 2]-na_wins[win_num, 1]+1) > int.win.max) next
			
					# Find window of values around NA window
					fit_win <- c(na_wins[win_num, 1] - fit.win.max, na_wins[win_num, 2] + fit.win.max)
				
					# Check that fit window doesn't exceed data boundaries
					if(fit_win[1] < 1) fit_win[1] <- 1
					if(fit_win[2] > dim(xyz)[3]) fit_win[2] <- dim(xyz)[3]
				
					# Check that window doesn't overlap into other NA values
					if(win_num > 1){
					#	if(fit_win[1] <= na_wins[win_num-1, 2]) fit_win[1] <- na_wins[win_num-1, 2] + 1
					}
					if(win_num < nrow(na_wins)){
					#	if(fit_win[2] >= na_wins[win_num+1, 1]) fit_win[2] <- na_wins[win_num+1, 1] - 1
					}
				
					# Set fit time points pre and post
					fit_itr_pre <- fit_win[1]:(na_wins[win_num, 1]-1)
					fit_itr_pos <- (na_wins[win_num, 2]+1):fit_win[2]
				
					# Skip if either is below min
					if(length(fit_itr_pre) < fit.win.min) next
					if(length(fit_itr_pos) < fit.win.min) next
				
					# Set that point has some interpolation
					point_has_inter[i] <- TRUE
				
					# Set iterations for full fit window
					full_fit_itr <- fit_itr_pre[1]:tail(fit_itr_pos,1)
				
					# Set indices that are NA
					na_inter_ind <- (length(fit_itr_pre)+1):(length(full_fit_itr)-length(fit_itr_pos))

					# Set span
					span <- round((tail(span.factor,1)/length(full_fit_itr)) + 0.01, 3)

					# For each dimension
					for(j in 1:dim(xyz)[2]){

						# Create dataframe
						data_frame <- data.frame(y=xyz[i, j, full_fit_itr], x=full_fit_itr)
					
						# Fit loess (higher span value corresponds to lower parameter fit)
						loess_main <- suppressWarnings(loess(y ~ x, data=data_frame, span=span))
					
						# Create smoothed points
						smoothed <- predict(loess_main, data_frame)
					
						# Save which indices were interpolated
						if(j == 1) inter_indices[[i]][[win_num]] <- na_wins[win_num, 1]:na_wins[win_num, 2]
					
						# Fill in smoothed values
						xyz_inter[i,j,na_wins[win_num, 1]:na_wins[win_num, 2]] <- smoothed[na_inter_ind]
					}
				}
			}
		}
		
		xyz_label <- c('x', 'y', 'z')
		dark_shade <- 130
		cols <- c(rgb(1,0,0), rgb(dark_shade/255,0,0), rgb(0,1,0), rgb(0,dark_shade/255,0), rgb(0,1,1), rgb(0,dark_shade/255,dark_shade/255))
		cols_sd <- c("#984EA3","#FF7F00","#FFFF33")

		### Plot
		if(!is.null(plot.diag) && any(point_has_inter)){

			# Change to indices
			point_has_inter <- which(point_has_inter)

			plot_height <- 3.8*length(point_has_inter)
			plot_width <- max(round(dim(xyz)[3]/150), 15)

			# Create layout matrix
			layout_mat <- cbind(1:length(point_has_inter))

			# If plot filename doesn't end in pdf, add
			if(!grepl('[.]pdf$', plot.diag, ignore.case=TRUE)) plot.diag <- paste0(plot.diag, '.pdf')

			pdf(plot.diag, height=plot_height, width=plot_width)
			layout(layout_mat)
			par(mar=c(4.5,4.5,3,1))
			xlim <- range(times)

			# For each point
			for(i in point_has_inter){
				
				# Re-order from start to end
				start_idx <- rep(NA, length(inter_indices[[i]]))
				for(k in 1:length(inter_indices[[i]])){
					start_idx[k] <- inter_indices[[i]][[k]][1]
				}
				inter_indices[[i]] <- inter_indices[[i]][order(start_idx)]

				# Get mean centered coordinates
				xyz_mc <- t(xyz_inter[i,,])
				for(j in 1:dim(xyz)[2]){
					xyz_mc[,j] <- xyz_inter[i,j,]-mean(xyz_inter[i,j,], na.rm=TRUE)
				}

				# Create plot window
				ylim <- range(xyz_mc, na.rm=TRUE)

				main_text <- paste0(dimnames(xyz)[[1]][i], ' (span: ', span, ')')

				plot(x=xlim, y=ylim, type='n', main=main_text, xlab='Time', ylab='Mean-centered values')

				for(j in 1:dim(xyz)[2]){

					for(k in 1:length(inter_indices[[i]])){
						
						# Skip if NULL
						if(is.null(inter_indices[[i]][[k]])) next
						
						if(fit_single_poly[i]){

							# Plot smoothed mean centered coordinates
							points(x=times[inter_indices[[i]][[k]]], y=xyz_mc[inter_indices[[i]][[k]],j], type='l', col=cols[j*2-1])
							points(x=times, y=xyz[i, j, ] - mean(xyz_inter[i,j,], na.rm=TRUE), col=cols[j*2], cex=0.8)
							points(x=times[inter_indices[[i]][[k]]], y=xyz_mc[inter_indices[[i]][[k]],j], pch=16, col=cols[j*2-1], cex=0.5)

						}else{

							if(k == 1){
								pre_inter <- 1:(inter_indices[[i]][[k]][1] - 1)
							}else{
								pre_inter <- (last_inter+1):(inter_indices[[i]][[k]][1] - 1)
							}
						
							# Plot smoothed mean centered coordinates
							points(x=times[pre_inter], xyz[i,j,pre_inter]-mean(xyz_inter[i,j,], na.rm=TRUE), type='l', col=cols[j*2])

							# Plot mean centered coordinates with interpolation
							#print(inter_indices[[i]][[k]])
							if(length(inter_indices[[i]][[k]]) == 1){
								points(x=times[inter_indices[[i]][[k]]], y=xyz_mc[inter_indices[[i]][[k]],j], pch=16, col=cols[j*2-1], cex=0.5)
							}else{
								points(x=times[inter_indices[[i]][[k]]], y=xyz_mc[inter_indices[[i]][[k]],j], type='l', lwd=2, col=cols[j*2-1], cex=0.7)
							}
						
							# Set last interpolated index
							last_inter <- tail(inter_indices[[i]][[k]], 1)
						}
					}

					if(!fit_single_poly[i]){
						# Plot anything after last interpolation
						pre_inter <- (last_inter+1):dim(xyz)[3]
						points(x=times[pre_inter], xyz[i,j,pre_inter]-mean(xyz_inter[i,j,], na.rm=TRUE), type='l', col=cols[j*2])
					}
				}
			}
			
			dev.off()
		}
	}
	
	## Interpolate transformation array?

	if(input_coordinates) return(xyz_inter)
	
	motion$xyz <- xyz_inter

	return(motion)
}