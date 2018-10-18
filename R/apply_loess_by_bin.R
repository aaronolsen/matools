apply_loess_by_bin <- function(x, bins, span, overlap=0){
	
	xf <- rep(NA, length(x))
	
	fitmat <- matrix(NA, nrow=0, ncol=length(x))
	fitbin <- c()

	i <- 1
	while(i < length(x)){
		
		# Get bin size
		first_not <- i + which(bins[i] != bins[i:length(bins)])[1] - 1

		# Reached the last segment, break
		if(is.na(first_not)) first_not <- length(bins) + 1
		
		x_idx <- max(1, i-overlap):min(first_not-1+overlap, length(x))
		x_idx_no <- i:(first_not-1)
		
		y <- x[x_idx]
		
		if(span[bins[i]] == 0){

			xf[x_idx_no] <- y[x_idx_no]
			fitmat <- rbind(fitmat, rep(NA, length(x)))

			fitmat[nrow(fitmat), x_idx] <- y
			fitbin <- c(fitbin, bins[i])

		}else{

			lp_loess <- suppressWarnings(loess(y ~ t, data=data.frame(t=x_idx, y=y), 
				span=(span[bins[i]]/length(y))+0.01))
				#span=0.07))
		
			xf[x_idx_no] <- predict(lp_loess, x_idx_no)
			fitmat <- rbind(fitmat, rep(NA, length(x)))

			fitmat[nrow(fitmat), x_idx] <- predict(lp_loess, x_idx)
			fitbin <- c(fitbin, bins[i])
		}

		i <- first_not
	}
	
	# Find overlap points
	overlap_cols <- which(colSums(!is.na(fitmat)) > 1)

	xf <- colMeans(fitmat, na.rm=TRUE)

	# Smooth overlap points
	if(length(overlap_cols) > 0){
		for(i in 1:length(overlap_cols)){

			if(overlap_cols[i] %in% c(1,length(xf))) next

			for(j in c(-1,1,0)){
				xf[overlap_cols[i]+j] <- mean(xf[overlap_cols[i]+c(-1,0,1)+j])
			}
		}
	}

	list(
		'x'=xf,
		'fitmat'=fitmat,
		'fitbin'=fitbin,
		'overlap'=overlap_cols
	)
}