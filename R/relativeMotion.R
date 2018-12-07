relativeMotion <- function(motion, fixed, ref.iter = 1){

	# fixed can be:
	#	names of xyz points
	#	future: xyz coordinates (matrix)
	#	future: body name (using tmat)

	if(!is.null(motion$xyz)){
	
		# Get xyz
		xyz <- motion$xyz

		# Set fixed coordinate reference
		if(is.vector(fixed)) ref_xyz <- xyz[fixed, , ref.iter]
		if(is.matrix(fixed)) ref_xyz <- fixed
		
		# Check for NA values
		if(sum(!is.na(ref_xyz[, 1])) < 3) stop('Fixed points must have at least 3 non-NA values at reference iteration. Input points have ', sum(!is.na(ref_xyz[, 1])), ' non-NA value(s) at reference iteration.')

		# Get number of iterations
		n_iter <- motion$n.iter

		# Each iteration
		for(iter in 1:n_iter){

			# Find translation and rotation to align coor to pts at point 1
			best_align <- bestAlign(ref_xyz, xyz[rownames(ref_xyz), , iter], xyz[, , iter])

			# Save copied alignment
			xyz[, , iter] <- best_align$mc

			# Apply transformation to transformation matrices if present
			if(!is.null(motion$tmat)) motion$tmat[, , , iter] <- applyTransform(motion$tmat[, , , iter], best_align$tmat)
		}
		
		# Replace xyz
		motion$xyz <- xyz
	}
	
	motion
}