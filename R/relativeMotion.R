relativeMotion <- function(motion, fixed, ref.iter = 1){

	# fixed can be:
	#	names of xyz points
	#	future: xyz coordinates (matrix)
	#	future: body name (using tmat)

	# Get number of iterations
	n_iter <- motion$n.iter
	
	# Check if there are transformations and xyz
	has_tmat <- ifelse(is.null(motion$tmat), FALSE, TRUE)
	has_xyz <- ifelse(is.null(motion$xyz), FALSE, TRUE)

	# Get input type
	if(is.vector(fixed) && length(fixed) == 1){
		fixed_type <- 'body_name'
	}else if(is.vector(fixed) && length(fixed) > 1){
		fixed_type <- 'xyz_names'
	}else if(is.matrix(fixed)){
		if(ncol(fixed) == 3){
			fixed_type <- 'xyz'
		}else if(ncol(fixed) == 4){
			fixed_type <- 'tmat_mat'
		}
	}else{
		fixed_type <- 'tmat_arr'
	}
	
	## Check inputs
	if(fixed_type == 'body_name'){

		# Check tmat for body name input
		if(!has_tmat) stop('If "fixed" input parameter is a body name then "motion" input must have transformations ("tmat").')

		# Check that body has transformations
		if(!fixed %in% dimnames(motion$tmat)[[3]]) stop('"Fixed" input , "', fixed, '", not found in motion$tmat.')
	}
	
	# Check xyz input
	if(fixed_type %in% c('xyz', 'xyz_names') && !has_xyz) stop('If "fixed" input parameter are xyz coordinates (or names of xyz coordinates) then "motion" must have an xyz object.')

	if(fixed_type %in% c('body_name', 'tmat_mat', 'tmat_arr')){
		
		# Apply transformation at each iteration
		for(iter in 1:n_iter){
		
			# Get inverse transformation of fixed body
			if(fixed_type == 'body_name'){
				inv_tmat <- solve(motion$tmat[, , fixed, iter])
			}else if(fixed_type == 'tmat_arr'){
				inv_tmat <- solve(fixed[, , iter])
			}else{
				inv_tmat <- solve(fixed)
			}

			# Apply transformation to transformation matrices if present
			if(has_tmat) motion$tmat[, , , iter] <- applyTransform(motion$tmat[, , , iter], inv_tmat)
			
			# Transform xyz coordinates
			if(has_xyz) motion$xyz[, , iter] <- applyTransform(motion$xyz[, , iter], inv_tmat)
		}

	}else{

		## Set relative motion using coordinates
		# Get xyz
		xyz <- motion$xyz

		# Set fixed coordinate reference
		if(is.vector(fixed) && length(fixed) > 1) ref_xyz <- xyz[fixed, , ref.iter]
		if(is.matrix(fixed)) ref_xyz <- fixed

		# Check for NA values
		if(sum(!is.na(ref_xyz[, 1])) < 3) stop('Fixed points must have at least 3 non-NA values at reference iteration. Input points have ', sum(!is.na(ref_xyz[, 1])), ' non-NA value(s) at reference iteration.')

		# Each iteration
		for(iter in 1:n_iter){

			# Find translation and rotation to align coor to pts at point 1
			best_align <- bestAlign(ref_xyz, xyz[rownames(ref_xyz), , iter], xyz[, , iter])

			# Save copied alignment
			xyz[, , iter] <- best_align$mc

			# Apply transformation to transformation matrices if present
			if(has_tmat) motion$tmat[, , , iter] <- applyTransform(motion$tmat[, , , iter], best_align$tmat)
		}
	
		# Replace xyz
		motion$xyz <- xyz
	}
	
	motion
}