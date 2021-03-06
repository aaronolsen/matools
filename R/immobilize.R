immobilize <- function(coor, ref){

	## Fixes motion of coor over multiple iterations relative to ref points (single iteration)
	
	# Get number of iterations
	n_iter <- dim(coor)[3]
	
	# Each iteration
	for(iter in 1:n_iter){

		# Find translation and rotation to align coor to pts at point 1
		coor[, , iter] <- bestAlign(ref, coor[rownames(ref), , iter], coor[, , iter])$mc
	}

	coor
}