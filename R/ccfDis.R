ccfDis <- function(x, y, lag.max, sets = rep(1, length(x)), correlation = c('ccf', 'pearson')){

	# Set lags
	lags <- -lag.max:lag.max
	
	# Get number of sets
	unique_sets <- unique(sets)

	# Set buffered lengths
	buff_len <- length(x) + lag.max*(length(unique_sets)-1)
	
	# Create buffered vectors
	x_buff <- rep(NA, buff_len)
	y_buff <- rep(NA, buff_len)
	
	# Fill buffered vectors
	n <- 1
	for(i in 1:length(unique_sets)){

		# Get indices in set
		which_in_set <- sets == unique_sets[i]
		
		# Get where to add to buffered vector
		add_at <- n:(n-1+sum(which_in_set))

		# Fill values
		x_buff[add_at] <- x[which_in_set]
		y_buff[add_at] <- y[which_in_set]
		
		# Update start index
		n <- tail(add_at, 1) + lag.max + 1
	}

	#print(rbind(x_buff, y_buff))
	
	acf_vals <- rep(NA, length(lags))
	
	for(i in 1:length(lags)){

		x_slide <- y_slide <- rep(NA, buff_len + abs(lags[i]))

		if(lags[i] < 0){
			x_slide[(abs(lags[i])+1):(buff_len + abs(lags[i]))] <- x_buff
			y_slide[1:buff_len] <- y_buff
		}else{
			x_slide[1:buff_len] <- x_buff
			y_slide[(abs(lags[i])+1):(buff_len + abs(lags[i]))] <- y_buff
		}

		# Calculate cross-correlation from formula
		# 	This matches ccf function with single continuous time series
		x_cc <- x_slide - mean(x_slide, na.rm=TRUE)
		y_cc <- y_slide - mean(y_slide, na.rm=TRUE)

		if(correlation[1] == 'ccf'){
			acf_vals[i] <- sum(x_cc * y_cc, na.rm=TRUE) / (sqrt(sum(x_cc^2, na.rm=TRUE)) * sqrt(sum(y_cc^2, na.rm=TRUE)))
		}else{
			acf_vals[i] <- cor(x_slide, y_slide, use='na.or.complete')
		}
	}
	
	# Find the lag at which CC is greatest
	which_abs_max <- which.max(abs(acf_vals))
	
	list('acf'=acf_vals, 'lag'=lags, 'acf.abs.max'=acf_vals[which_abs_max], 'lag.abs.max'=lags[which_abs_max])
}