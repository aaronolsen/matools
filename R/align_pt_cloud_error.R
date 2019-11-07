align_pt_cloud_error <- function(p, m1, m2, center=c(0,0,0)){

	# Set transformation
	tmat1 <- tmat2 <- tmat3 <- diag(4)
	tmat1[1:3, 4] <- center+p[1:3]
	tmat2[1:3, 1:3] <- rotationMatrixZYX_ma(p[4:6])
	tmat3[1:3, 4] <- -center
	tmat <- tmat1 %*% tmat2 %*% tmat3

	# Transform m2 by tmat = m2_t
	m2_t <- mtransform(m2, tmat)

	# For each point in m2_t find closest point in m1
	if(FALSE){
		d_pts <- rep(NA, length=nrow(m2_t))
		for(i in 1:nrow(m2_t)){

			# Find mean square difference from m2 point to each point in m1
			pt_to_pt <- rep(NA, length=nrow(m1))
			#for(j in 1:nrow(m1)) pt_to_pt[j] <- mean((m2_t[i,]-m1[j,])^2)
			for(j in 1:nrow(m1)) pt_to_pt[j] <- sum(abs(m2_t[i,]-m1[j,]))

			# Save minimum
			d_pts[i] <- min(pt_to_pt)
		}

		return(sum(d_pts))
	}
	
	#
	min_pt_cloud_dist(m1, m2_t)
}