ref_rotate_error <- function(p, ref.points, fit.points, center, axis = NULL){

	# Create transformation matrix
	tmat1 <- tmat2 <- tmat3 <- diag(4)
	tmat1[1:3,4] <- center
	
	if(!is.null(axis)){
		tmat2[1:3,1:3] <- tMatrixEP_ma(axis, p)
	}else{
		tmat2[1:3,1:3] <- rotationMatrixZYX_ma(p)
	}

	tmat3[1:3,4] <- -center

	# Save transformation matrix
	tmat <- tmat1 %*% tmat2 %*% tmat3

	# Apply test transformation to reference points
	ref.points <- applyTransform(ref.points, tmat)	

	# Return error
	return(sqrt(mean((ref.points - fit.points)^2)))
}