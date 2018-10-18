fourierTransform <- function(x, y, bins = NULL, freq.norm.min = NULL, freq.norm.max = NULL){

	# 
	#bins <- bins + 2

	# Compute Fourier transform
	y_ft <- fft(y)

	# Get amplitudes
	amps <- 2*(Mod(y_ft) / length(y_ft))

	# Get frequencies
	freqs <- 0:(length(x)-1)

	# Get indices to keep
	if(length(amps) %% 2 == 1){	# odd
		half1_idx <- 2:floor((length(amps)/2))
		half2_idx <- length(amps):(floor((length(amps)/2))+3)
	}else{
		half1_idx <- 2:floor((length(amps)/2))
		half2_idx <- length(amps):(floor((length(amps)/2))+2)
	}

	# Remove duplicate values
	amps <- amps[half1_idx]
	freqs <- freqs[half1_idx]
	
	# Get normalized frequencies
	freqs_norm <- freqs / diff(range(x))

	# Apply any cut-offs
	if(!is.null(freq.norm.min)){
		ind_keep <- freqs_norm >= freq.norm.min
		freqs_norm <- freqs_norm[ind_keep]
		amps <- amps[ind_keep]
		freqs <- freqs[ind_keep]
		y_ft <- y_ft[ind_keep]
	}
	if(!is.null(freq.norm.max)){
		ind_keep <- freqs_norm <= freq.norm.max
		freqs_norm <- freqs_norm[ind_keep]
		amps <- amps[ind_keep]
		freqs <- freqs[ind_keep]
		y_ft <- y_ft[ind_keep]
	}
	
	if(!is.null(bins) && length(freqs_norm) < bins) warning(paste0("Number of frequencies (", length(freqs_norm), ") is less than the number of bins (", bins, ")."))

	if(!is.null(bins)){
	
		bin_limits <- seq(0, (length(amps)+1), length=bins+1)

		y_ft_binned <- rep(NA, bins)
		amps_binned <- rep(NA, bins)
		freqs_binned <- rep(NA, bins)
		freqs_norm_binned <- rep(NA, bins)

		max_last <- 0

		for(i in 1:(length(bin_limits)-1)){
		
			ind1 <- ceiling(bin_limits[i])
			ind2 <- floor(bin_limits[i+1])

			# Get indices for bin
			indices <- ind1:ind2
			indices <- indices[indices > 0]
			indices <- indices[indices <= length(amps)]
			
			# Make sure indices are all greater than previous maximum index
			indices <- indices[indices > max_last]

			# Get binned values
			y_ft_binned[i] <- sum(y_ft[indices])
			amps_binned[i] <- sum(amps[indices])
			freqs_binned[i] <- mean(freqs[indices])
			freqs_norm_binned[i] <- mean(freqs_norm[indices])
			
			max_last <- max(indices)
		}
		
		y_ft <- y_ft_binned
		amps <- amps_binned
		freqs <- freqs_binned
		freqs_norm <- freqs_norm_binned
	}

	#print(amps)
	#hist_amp <- hist(amps[half1_idx], breaks=10, plot=FALSE)
	#hist_amp_norm <- hist_amp$counts*hist_amp$mids
	#plot(hist_amp$mids, hist_amp_norm, type='h', cex=10)
	#print(half1_idx)
	
	#plot(freqs_norm_binned[half1_idx], amps[half1_idx], type='l')
	
	# Set baseline/constant
	constant <- Re(y_ft[1]) / length(y_ft)

	rlist <- list(
		#'ft' = y_ft,
		'freqs' = freqs,
		'freqs.norm' = freqs_norm,
		'amps' = amps
		#'constant' = constant,
		#'half1.idx' = half1_idx,
		#'half2.idx' = half2_idx
	)
	
	return(rlist)
}