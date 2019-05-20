ref_rotate_error <- function(p, ref.points, fit.points, center, axis = NULL, 
	ref.plane.point = NULL, fit.plane.point = NULL, fit.plane.normal = NULL, plane.marker.idx = NULL){

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

	# Add projection of last ref point into plane as last (corresponding) fit point
	if(!is.null(fit.plane.point)){
		if(is.matrix(ref.points)){
			fit.points <- rbind(fit.points, pointPlaneProj_ma(ref.points[plane.marker.idx, ], fit.plane.point, fit.plane.normal))
		}else{
			fit.points <- rbind(fit.points, pointPlaneProj_ma(ref.points, fit.plane.point, fit.plane.normal))
		}
	}

	# Return error
	return(sqrt(mean((ref.points - fit.points)^2)))
}