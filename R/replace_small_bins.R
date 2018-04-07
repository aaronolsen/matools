replace_small_bins <- function(x, min.size, replace.with='max'){

	i <- 1
	while(i < length(x)){
		
		# Get bin size
		first_not <- i + which(x[i] != x[i:length(x)])[1] - 1

		if(is.na(first_not)) first_not <- length(x) + 1

		bin_size <- first_not-i
		
		#
		if(bin_size < min.size){

			# Replace values
			x[i:(first_not-1)] <- do.call(replace.with, list(x[c(max(1, i-1), min(first_not, length(x)))]))
		}

		i <- first_not
	}
	
	x
}