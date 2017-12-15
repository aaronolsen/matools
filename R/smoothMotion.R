smoothMotion <- function(xyz, span.factor = 18, min.times = 10, plot.diag = NULL, diag.combine.dim = TRUE){

	xyz_label <- c('x', 'y', 'z')
	cols <- c(rgb(1,0,0), rgb(180/255,0,0), rgb(0,1,0), rgb(0,180/255,0), rgb(0,0,1), rgb(0,0,180/255))

	if(length(dim(xyz)) == 3){

		# Time vector
		times <- 1:dim(xyz)[3]
		
		# Create smoothed array
		xyz_smooth <- xyz
		
		if(!is.null(plot.diag)){
			if(diag.combine.dim){
				plot_height <- 2.8*dim(xyz)[1]
				plot_width <- max(round(dim(xyz)[3]/230), 9)
				layout_mat <- cbind(1:dim(xyz)[1])
			}else{
				plot_height <- 3*dim(xyz)[1]
				plot_width <- max(round(dim(xyz)[3]/70), 12)
				layout_mat <- matrix(1:(dim(xyz)[1]*dim(xyz)[2]), ncol=3, byrow=TRUE)
			}

			# If plot filename doesn't end in pdf, add
			if(!grepl('[.]pdf$', plot.diag, ignore.case=TRUE)) plot.diag <- paste0(plot.diag, '.pdf')

			pdf(plot.diag, height=plot_height, width=plot_width)
			layout(layout_mat)
			par(mar=c(4.5,4.5,3,1))
			xlim <- range(times)
		}
		
		# Smooth each point
		for(i in 1:dim(xyz)[1]){
		
			# Get number of non-NA frames
			num_non_NA <- sum(!is.na(xyz[i, 1, ]))

			if(num_non_NA < min.times){
				if(!is.null(plot.diag)){
					if(diag.combine.dim){
						plot(xlim, c(0,1), type='n', main=dimnames(xyz)[[1]][i], xlab='Time', ylab='Mean-centered values')
					}else{
						plot(xlim, c(0,1), type='n', main=dimnames(xyz)[[1]][i], xlab='Time', ylab='x-value')
						plot(xlim, c(0,1), type='n', main=dimnames(xyz)[[1]][i], xlab='Time', ylab='y-value')
						plot(xlim, c(0,1), type='n', main=dimnames(xyz)[[1]][i], xlab='Time', ylab='z-value')
					}
				}
				next
			}

			# Determine span based on number of non-NA frames
			# Higher span.factor = more smoothing
			span <- round((span.factor/num_non_NA) + 0.01, 3)

			if(!is.null(plot.diag)){
				if(diag.combine.dim){

					# Get mean centered coordinates
					xyz_mc <- t(xyz[i,,])
					for(j in 1:ncol(xyz_mc)){
						xyz_mc[,j] <- xyz[i,j,]-mean(xyz[i,j,], na.rm=TRUE)
					}

					# Create plot window
					plot(x=xlim, y=range(xyz_mc, na.rm=TRUE), type='n', main=paste0(dimnames(xyz)[[1]][i], 
						' (span: ', span, ')'), xlab='Time', ylab='Mean-centered values')
				}
			}

			# Smooth each dimension
			for(j in 1:dim(xyz)[2]){

				# Create data frame
				data_frame <- data.frame(y=xyz[i,j,!is.na(xyz[i,j,])], x=times[!is.na(xyz[i,j,])])

				# Smooth - higher span value corresponds to smoother fit
				loess_fit <- loess(y ~ x, data=data_frame, span=span)

				# Create smoothed points
				smoothed <- predict(loess_fit, data_frame)

				if(!is.null(plot.diag)){
					if(diag.combine.dim){
						points(x=times, y=xyz_mc[,j], col=cols[j*2-1], cex=0.7)
						points(x=times[!is.na(xyz[i,j,])], smoothed-mean(xyz[i,j,], na.rm=TRUE), type='l', col=cols[j*2])
					}else{
						plot(x=times, y=xyz[i,j,], main=paste0(dimnames(xyz)[[1]][i], '.', xyz_label[j], 
							' (span: ', span, ')'), xlab='Time', ylab=paste0(xyz_label[j], '-value'), col=cols[j*2])
						points(x=times[!is.na(xyz[i,j,])], smoothed, type='l', col=cols[j*2-1])
					}
				}

				# Replace unsmoothed with smoothed coordinates
				xyz_smooth[i,j,!is.na(xyz[i,j,])] <- smoothed
			}
		}

		if(!is.null(plot.diag)) dev.off()

	}else{
		
		# Transformation array
	}
	
	xyz_smooth
}