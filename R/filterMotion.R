filterMotion <- function(motion, min.value, min.length, xyz.use = NULL, span.factor = 18, filter = c('mean_abs_deriv'), 
	return = c('motion', 'remove.rows'), plot.diag = NULL){
	
	# Standardize filter capitalization
	filter <- tolower(filter)
	
	# Filter by mean absolute derivative
	if(filter == 'mean_abs_deriv'){
	
		# Check that xyz.use has values
		if(!is.null(xyz.use) && length(xyz.use) == 0) stop("length(xyz.use) is 0.")
		
		# Set points to use in filtering
		if(is.null(xyz.use)) xyz.use <- 1:dim(motion$xyz)[1]
		
		# Get points to filter
		xyz <- motion$xyz[xyz.use, , ]
		
		# Smooth fit points
		# Higher span.factor = more smoothing
		if(!is.null(span.factor)){
			xyz_smooth <- smoothMotion(xyz, span.factor=span.factor, plot.diag=NULL)
		}else{
			xyz_smooth <- xyz
		}

		# Create matrix for point derivative
		xyz_d <- matrix(NA, nrow=dim(xyz)[1], ncol=dim(xyz)[3], dimnames=list(dimnames(xyz)[[1]], NULL))

		# Create matrix for mean absolute derivative
		xyz_mad <- rep(NA, dim(xyz)[3])

		# Set times/iterations
		if('time' %in% names(motion)){
			times <- motion[['time']]
			time_label <- 'Time'
			
			# Convert min.length to iterations
			frame_rate <- dim(xyz)[3] / max(times)
			min.length <- min.length * frame_rate

		}else{
			times <- 1:dim(xyz)[3]
			time_label <- 'Iterations'
		}
		
		# Create vector of time-to-time distances
		dist_xyz <- rep(NA, dim(xyz)[3])
		diff_xyz <- rep(NA, dim(xyz)[3])

		# Find derivative of points
		for(i in 1:dim(xyz)[1]){

			#
			diff_xyz <- distPointToPoint(t(xyz_smooth[i,,]))
			
			# Find derivative
			xyz_d[i, ] <- colMeans(rbind(c(NA, diff_xyz/diff(times)),
				c(diff_xyz/diff(times), NA)), na.rm=TRUE)
		}
		
		# Find mean absolute derivative
		xyz_mad <- apply(abs(xyz_d), 2, 'max', na.rm=TRUE)
		
		# Find where values are under threshold
		min_mad_under <- xyz_mad < min.value
		
		# Check whether there are any values above threshold
		if(sum(min_mad_under) > 0){

			# Find starts and ends of TRUE
			true_start <- which(diff(c(FALSE, min_mad_under, FALSE)) == 1)
			true_end <- which(diff(c(FALSE, min_mad_under, FALSE)) == -1)

			for(i in 1:length(true_start)){
			
				# Find
				start_at <- true_start[i]
				end_at <- true_end[i]-1

				# If length is less than minimum, change to false
				if(((end_at+1)-start_at) < min.length) min_mad_under[start_at:end_at] <- FALSE
			}
		
			# Set indices to remove
			remove_row <- min_mad_under
			
			# Set trim
			trim <- TRUE

		}else{

			# Set trim
			trim <- FALSE
			remove_row <- rep(FALSE, dim(xyz)[3])
			xyz_trim <- xyz
		}

		# Make diagnostic plot
		if(!is.null(plot.diag)){
		
			# Get filtered points
			xyz_trim <- xyz[, , !remove_row]

			# Set colors
			cols <- rainbow(n=dim(xyz)[1], start=0, end=0.7)

			# If plot filename doesn't end in pdf, add
			if(!grepl('[.]pdf$', plot.diag, ignore.case=TRUE)) plot.diag <- paste0(plot.diag, '.pdf')
			
			# Find new starts and ends of TRUE
			true_start <- which(diff(c(FALSE, min_mad_under, FALSE)) == 1)
			true_end <- which(diff(c(FALSE, min_mad_under, FALSE)) == -1)
			
			# Set plot width, no greater than 15
			plot_width <- min(max(6, round(1600/300)), 15)

			pdf(plot.diag, width=plot_width, height=11)

			layout(cbind(1:4), height=c(1,1,1,1))
			par(mar=c(4.5,5.4,2,2))

			##
			# Plot point derivative
			xlim <- range(times)
			ylim <- c(0, max(xyz_d, na.rm=TRUE))
			plot(xlim, ylim, type='n', ylab='Derivative of point position', xlab=time_label)
			abline(h=0, lty=2, col=gray(0.7))
			for(i in 1:nrow(xyz_d)){
				points(times, xyz_d[i, ], type='l', col=cols[i])
			}

			##
			# Plot mean absolute derivative
			xlim <- range(times)
			ylim <- c(0, max(xyz_mad, na.rm=TRUE))
			plot(xlim, ylim, type='n', ylab='Mean absolute value of derivative', xlab=time_label)
			#polygon(x=rbind(c(xlim[1],min.value),c(xlim[2],min.value),c(xlim[2],0),c(xlim[1],0)), 
			#	col=rgb(0.8,0.8,1), border=NA)

			# Plot polygon in background for removed portions
			if(length(true_start) > 0){
				for(i in 1:length(true_start)){
			
					# Find
					start_at <- times[true_start[i]]
					end_at <- times[true_end[i]-1]

					polygon(x=rbind(c(start_at,ylim[2]),c(end_at,ylim[2]),c(end_at,ylim[1]),c(start_at,ylim[1])), 
						col=rgb(0.9,0.9,0.9), border=NA)
				}
			}

			abline(h=0, lty=2, col=gray(0.7))
			points(times, xyz_mad, type='l')

			text(x=xlim[2], y=min.value+0.03*diff(ylim), labels='min.value', pos=2, col=gray(0.7))
			abline(h=min.value, lty=3, col=gray(0.7))

			##
			# Get mean-centered position for each point
			xyz_mc <- xyz
			for(i in 1:dim(xyz)[1]) for(j in 1:dim(xyz)[2]) xyz_mc[i, j, ] <- xyz[i, j, ]-mean(xyz[i, j, ], na.rm=TRUE)

			# Get smoothed mean-centered position for each point
			xyz_s_mc <- xyz_smooth
			for(i in 1:dim(xyz)[1]) for(j in 1:dim(xyz)[2]) xyz_s_mc[i, j, ] <- xyz_smooth[i, j, ]-mean(xyz_smooth[i, j, ], na.rm=TRUE)
			
			# Plot point positions over times
			xlim <- range(times)
			ylim <- range(xyz_mc, na.rm=TRUE)
			plot(xlim, ylim, type='n', ylab='Mean-centered position', xlab=time_label)

			# Plot polygon in background for removed portions
			if(length(true_start) > 0){
				for(i in 1:length(true_start)){
			
					# Find
					start_at <- times[true_start[i]]
					end_at <- times[true_end[i]-1]

					polygon(x=rbind(c(start_at,ylim[2]),c(end_at,ylim[2]),c(end_at,ylim[1]),c(start_at,ylim[1])), 
						col=rgb(0.9,0.9,0.9), border=NA)
				}
			}
		
			# Plot mean-centered points
			for(i in 1:dim(xyz_mc)[1]) for(j in 1:dim(xyz_mc)[2]) points(times, xyz_mc[i, j, ], cex=0.2, col=cols[i])

			# Plot smoothed mean-centered points
			for(i in 1:dim(xyz_s_mc)[1]) for(j in 1:dim(xyz_s_mc)[2]) points(times, xyz_s_mc[i, j, ], type='l', lwd=1, col=cols[i])
		

			##
			# Get mean-centered position for each point
			xyz_trim_mc <- xyz_trim
			for(i in 1:dim(xyz_trim)[1]) for(j in 1:dim(xyz_trim)[2]) xyz_trim_mc[i, j, ] <- xyz_trim[i, j, ]-mean(xyz_trim[i, j, ], na.rm=TRUE)

			# Plot point positions over times
			xlim <- c(1,dim(xyz_trim_mc)[3])
			ylim <- range(xyz_trim_mc, na.rm=TRUE)
			plot(xlim, ylim, type='n', ylab='Mean-centered position (trimmed)', xlab=time_label)
		
			for(i in 1:dim(xyz_trim_mc)[1]) for(j in 1:dim(xyz_trim_mc)[2]) points(xlim[1]:xlim[2], xyz_trim_mc[i, j, ], type='l', col=cols[i])

			dev.off()
		}
	}
	
	# Remove rows
	remove_rows <- which(remove_row)
	if(return[1] == 'motion') return(removeRows(motion, remove.rows=remove_rows))
	if(return[1] == 'remove.rows') return(remove_rows)
}