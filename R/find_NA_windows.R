find_NA_windows <- function(x, sort = FALSE){

	# Skip if all values are not NA
	if(!any(is.na(x))) return(NULL)

	# Skip if all values are NA
	if(!any(!is.na(x))) return(NULL)

	# Get which values are NA
	which_na <- which(is.na(x))
	which_nna <- which(!is.na(x))

	# Get first na
	first_na <- which_na[1]
	
	# Sequence starts with NA
	if(first_na == 1){
		
		# Find first NA after non NA
		first_na <- which_na[which_na > which_nna[1]][1]
		
		# If none, return NULL
		if(is.na(first_na)) return(NULL)
	}

	# Find first non NA after first NA
	last_na <- which_nna[which_nna > first_na][1] - 1

	# If none, return NULL
	if(is.na(last_na)) return(NULL)

	# Create matrix
	rmat <- matrix(NA, nrow=1, ncol=2)
	
	# Add first instance to list
	rmat[1,] <- c(first_na, last_na)
	
	# Continue finding more
	at_end <- FALSE
	while(!at_end){
		
		# Find first NA after last na
		next_first_na <- which_na[which_na > last_na][1]
		
		# If none, return NULL
		if(is.na(next_first_na)) break
		
		# Find next last NA
		next_last_na <- which_nna[which_nna > next_first_na][1] - 1
		
		# If none, return NULL
		if(is.na(next_last_na)) break
		
		# Add next to list
		rmat <- rbind(rmat, c(next_first_na, next_last_na))
		
		# Update
		last_na <- next_last_na
	}
	
	if(sort){
		row_len <- rmat[,2]-rmat[,1]+1
		rmat <- rmat[order(row_len), , drop=FALSE]
	}
	
	rmat
}