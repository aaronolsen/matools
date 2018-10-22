smoothMotion <- function(motion, span.factor = 18, min.times = 10, plot.diag = NULL, n.bins = NULL, 
	adaptive = FALSE, smooth.bins = NULL, bin.replace.min = 10, sd.win = 10){
	# , span.bin.dec = 0.02){

	# Reverse span factor to coincide with smooth bins order
	span.factor <- rev(span.factor)

	# Get xyz coordinates
	input_coordinates <- FALSE
	if(is.list(motion)){
		xyz <- motion$xyz
	}else{
		input_coordinates <- TRUE
		xyz <- motion
	}
	
	xyz_label <- c('x', 'y', 'z')
	dark_shade <- 130
	cols <- c(rgb(1,0,0), rgb(dark_shade/255,0,0), rgb(0,1,0), rgb(0,dark_shade/255,0), rgb(0,1,1), rgb(0,dark_shade/255,dark_shade/255))
	cols_sd <- c("#984EA3","#FF7F00","#FFFF33")

	if(length(dim(xyz)) == 3){

		### Perform initial smoothing
		# Create arrays
		xyz_smooth <- xyz*NA
		xyz_dev <- xyz*NA
		xyz_dist <- xyz[,1,]*NA
		
		# Time vector
		times <- 1:dim(xyz)[3]

		# For each point
		for(i in 1:dim(xyz)[1]){

			# Get number of non-NA frames
			num_non_NA <- sum(!is.na(xyz[i, 1, ]))

			if(num_non_NA < min.times) next

			# Determine span based on number of non-NA frames
			# Higher span.factor = more smoothing
			span_main <- round((tail(span.factor,1)/num_non_NA) + 0.01, 3)

			# For each dimension
			for(j in 1:dim(xyz)[2]){
			
				# Create data frame
				data_frame <- data.frame(y=xyz[i,j,!is.na(xyz[i,j,])], x=times[!is.na(xyz[i,j,])])

				# Smooth - higher span value corresponds to smoother fit
				loess_main <- suppressWarnings(loess(y ~ x, data=data_frame, span=span_main))

				# Create smoothed points
				smoothed <- predict(loess_main, data_frame)
				
				# Find deviation from smoothed
				xyz_dev[i,j,!is.na(xyz[i,j,])] <- abs(smoothed-xyz[i,j,!is.na(xyz[i,j,])])
				
				# Replace unsmoothed with smoothed coordinates
				xyz_smooth[i,j,!is.na(xyz[i,j,])] <- smoothed
			}

			# Find distance from smoothed
			xyz_dist[i,!is.na(xyz[i,1,])] <- dppt(t(xyz_smooth[i,,!is.na(xyz[i,1,])]), t(xyz[i,,!is.na(xyz[i,1,])]))
		}
		
		### Apply adaptive smoothing
		if(adaptive){

			## Create smooth bins
			# Get range of deviations
			dev_range <- range(xyz_dev, na.rm=TRUE)
			dev_mean <- mean(xyz_dev, na.rm=TRUE)
			
			# Set number of bins from span.factor length if NULL
			if(is.null(n.bins)) n.bins <- length(span.factor)
			
			# Create smooth bins
			if(!is.null(smooth.bins)){
				smooth_bins <- smooth.bins
			}else{

				smooth_bins <- c(0, exp(1)^seq(log(dev_mean), log(dev_range[2]*1.01), length=n.bins))
				#print(smooth_bins)
				#smooth_bins <- c(0, seq(dev_mean, dev_range[2]*1.01, length=n.bins))
			}
			
			#print(smooth_bins)

			# Create bin matrix
			smooth_bin_mat <- matrix(NA, nrow=length(smooth_bins)-1, ncol=3)
			for(i in 1:(length(smooth_bins)-1)) smooth_bin_mat[i,] <- c(smooth_bins[i:(i+1)], i)

			# Set span factor for each bin
			if(length(span.factor) == 2){
				span_factors <- sort(round(seq(span.factor[1], span.factor[2], length=n.bins), 1), decreasing=TRUE)
			}else if(length(span.factor) == n.bins){
				span_factors <- sort(span.factor, decreasing=TRUE)
			}else{
				stop(paste0('If length of span.factor (', length(span.factor), ') is not 2 it must be equal to the number of bins (', n.bins, ')'))
			}

			xyz_dev_wins <- xyz*NA
			dev_wins_bin <- xyz*NA
			#overlap_pts <- list()

			# For each point
			for(i in 1:dim(xyz)[1]){

				# Find first and last non-na
				is_na <- which(!is.na(xyz[i,1,]))

				# All NA
				if(length(is_na) == 0) next

				# Get first and last
				first_last <- is_na[c(1, length(is_na))]
				
				# Set range
				frange <- first_last[1]:first_last[2]
				
				# Number of non-NA values is less than min.times
				if(length(frange) < min.times) next
				
				# Check for NAs within
				if(sum(is.na(xyz[i,1,frange]) > 0)) next

				# For each dimension
				for(j in 1:dim(xyz)[2]){

					# Find mean deviation by sliding window
					xyz_dev_wins[i,j,frange] <- abs(slide_fun(xyz_dev[i,j,frange], window=sd.win, step=1, fun='mean'))
				
					# Bin values
					dev_wins_bin[i,j,frange] <- bin_replace(xyz_dev_wins[i,j,frange], smooth_bin_mat)

					# Replace small bins with maximal value in range around bin
					dev_wins_bin[i,j,frange] <- replace_small_bins(dev_wins_bin[i,j,frange], min.size=bin.replace.min, replace.with='max')

					# Apply loess filter by bin
					apply_loess <- apply_loess_by_bin(xyz[i,j,frange], dev_wins_bin[i,j,frange], span_factors, overlap=1)
					xyz_smooth[i,j,frange] <- apply_loess$x
					#overlap_pts[[(i-1)*dim(xyz)[2]+j]] <- apply_loess$overlap
				}

				# Find distance from smoothed
				xyz_dist[i,frange] <- dppt(t(xyz_smooth[i,,frange]), t(xyz[i,,frange]))
			}
		}

		### Plot
		if(!is.null(plot.diag)){

			plot_height <- 3.6*dim(xyz)[1]
			layout_mat <- cbind(1:dim(xyz)[1])

			# Create layout matrix
			if(adaptive){
				layout_mat <- create_layout_matrix(num.elem=dim(xyz)[1], elem.plots=3, ncol=3)
				plot_width <- max(round(dim(xyz)[3]/130), 22)
			}else{
				layout_mat <- create_layout_matrix(num.elem=dim(xyz)[1], elem.plots=2, ncol=2)
				plot_width <- max(round(dim(xyz)[3]/150), 18)
			}

			# If plot filename doesn't end in pdf, add
			if(!grepl('[.]pdf$', plot.diag, ignore.case=TRUE)) plot.diag <- paste0(plot.diag, '.pdf')

			pdf(plot.diag, height=plot_height, width=plot_width)
			layout(layout_mat)
			par(mar=c(4.5,4.5,3,1))
			xlim <- range(times)

			# For each point
			for(i in 1:dim(xyz)[1]){

				#
				if(sum(!is.na(xyz_smooth[i,,])) == 0){
					if(!adaptive){
						plot(xlim, c(0,1), type='n', main=dimnames(xyz)[[1]][i], xlab='Time', ylab='Mean-centered values')
						plot(xlim, c(0,1), type='n', main=dimnames(xyz)[[1]][i], xlab='Time', ylab='Deviation')
					}else{
						plot(xlim, c(0,1), type='n', main=dimnames(xyz)[[1]][i], xlab='Time', ylab='Mean-centered values')
						plot(xlim, c(0,1), type='n', main=dimnames(xyz)[[1]][i], xlab='Time', ylab='Deviation by window')
						plot(xlim, c(0,1), type='n', main=dimnames(xyz)[[1]][i], xlab='Time', ylab='Deviation')
					}
					next
				}

				# Get mean centered coordinates
				xyz_mc <- t(xyz[i,,])
				for(j in 1:dim(xyz)[2]){
					xyz_mc[,j] <- xyz[i,j,]-mean(xyz[i,j,], na.rm=TRUE)
				}

				# Create plot window
				ylim <- range(xyz_mc, na.rm=TRUE)
				
				if(adaptive){
					# Find bins used
					bins_used <- sort(unique(c(dev_wins_bin[i,,])), decreasing=TRUE)
					if(length(bins_used) == 0) bins_used <- 1
					main_text <- paste0(dimnames(xyz)[[1]][i], ' (span factors used: ', paste0(span_factors[bins_used], collapse=','), ')')
				}else{
					main_text <- paste0(dimnames(xyz)[[1]][i], ' (span: ', span_main, ')')
				}
				
				plot(x=xlim, y=ylim, type='n', main=main_text, xlab='Time', ylab='Mean-centered values')

				for(j in 1:dim(xyz)[2]){

					if(adaptive){
						#if(length(overlap_pts[[(i-1)*dim(xyz)[2]+j]]) > 0){
						#	for(k in 1:length(overlap_pts[[(i-1)*dim(xyz)[2]+j]])){
						#		abline(v=overlap_pts[[(i-1)*dim(xyz)[2]+j]][k], lty=2, col=gray(0.8))
						#	}
						#}
					}

					# Plot raw mean centered coordinates
					points(x=times, y=xyz_mc[,j], type='l', lwd=2, col=cols[j*2-1], cex=0.7)

					# Plot smoothed mean centered coordinates
					points(x=times[!is.na(xyz[i,j,])], xyz_smooth[i,j,!is.na(xyz[i,j,])]-mean(xyz[i,j,], na.rm=TRUE), type='l', col=cols[j*2])
				}

				if(adaptive){

					if(sum(!is.na(xyz_dev_wins[i,,])) == 0){
						plot(xlim, c(0,1), type='n', main=dimnames(xyz)[[1]][i], xlab='Time', ylab='Initial smooth deviation by window')
					}else{

						main_text <- paste0(dimnames(xyz)[[1]][i], ' (mean: ', round(mean(xyz_dev_wins[i,,], na.rm=TRUE), 2), 
							'; max: ', round(max(xyz_dev_wins[i,,], na.rm=TRUE), 2), ')')

						plot_ylim <- range(xyz_dev_wins[i,,], na.rm=TRUE)

						plot(x=xlim, y=plot_ylim, type='n', main=main_text, xlab='Time', ylab='Initial smooth deviation by window')

						#				
						for(j in 1:dim(xyz)[2]){

							points(x=times, y=xyz_dev_wins[i,j,], type='l', col=cols[j*2-1])

							dev_wins_bin_norm <- (dev_wins_bin[i,j,] - 1) / (length(smooth_bins)-2)
							dev_wins_bin_norm_scale <- dev_wins_bin_norm*diff(plot_ylim) + plot_ylim[1]
							points(x=times, y=dev_wins_bin_norm_scale, type='l', lwd=0.5, col=cols[j*2])
						}
					}
				}

				# Plot distance from smoothed to original coordinates
				main_text <- paste0(dimnames(xyz)[[1]][i], ' (mean: ', round(mean(xyz_dist[i,], na.rm=TRUE), 2), 
					'; max: ', round(max(xyz_dist[i,], na.rm=TRUE), 2), ')')

				plot(x=xlim, y=range(xyz_dist[i,], na.rm=TRUE), type='n', main=main_text, xlab='Time', ylab='Deviation')

				points(x=times, y=xyz_dist[i,], type='l')
			}
		}

		if(!is.null(plot.diag)) dev.off()

	}else{
		
		# Transformation array
	}

	if(input_coordinates) return(xyz_smooth)
	
	motion$xyz <- xyz_smooth

	return(motion)
}