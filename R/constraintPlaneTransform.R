constraintPlaneTransform <- function(xyz, ids, child.type, axis.with.vp = TRUE){

	# axis.with.vp : whether to include virtual points (with beads) to fit axis

	dptp_mat <- matrix(NA, nrow=length(ids), ncol=7, dimnames=list(NULL, c('pre.min', 'pre.mean', 'pre.min', 'angle', 'post.min', 'post.mean', 'post.max')))

	# Set identity matrix as default
	tmat <- diag(4)

	i <- 1
	for(id in ids){

		# Get constraint plane vectors and points
		cpv <- sort(rownames(xyz)[grepl(paste0('_cp', id, 'v'), rownames(xyz))])
		cpp <- sort(rownames(xyz)[grepl(paste0('_cp', id, 'p'), rownames(xyz))])

		# Get bodies
		cpv_body <- strsplit(cpv, '_')[[1]][1]
		cpp_body <- strsplit(cpp, '_')[[1]][1]

		if(child.type == 'p'){

			# Get beads
			cpp_beads <- rownames(xyz)[grepl(paste0(cpp_body, '_bead'), rownames(xyz))]
			
			# Add virtual points
			if(axis.with.vp) cpp_beads <- c(cpp_beads, rownames(xyz)[grepl(paste0('-', cpp_body), rownames(xyz))])
			
			# Get xyz points associated with body
			cpp_body_xyz <- rownames(xyz)[grepl(paste0(cpp_body, '(_|-)'), rownames(xyz))]

			# Get plane normal vector
			plane_n <- uvector_ma(cprod(xyz[cpv[2],]-xyz[cpv[1],], xyz[cpv[3],]-xyz[cpv[1],]))

			# Find distance from points to plane (positive means on the side of the normal vector)
			dptp <- distPointToPlane(xyz[cpp,], n=plane_n, q=xyz[cpv[1],])
			dptp_mat[i, 1:3] <- c(min(dptp, na.rm=TRUE), mean(dptp, na.rm=TRUE), max(dptp, na.rm=TRUE))

			# If all values <= 0 skip
			if(sum(dptp < 0) == 0){
				i <- i + 1
				next
			}

			# Fit line to beads to find axis of rotation
			tmat_aor <- fitLine3D(xyz[cpp_beads, ])
			tmat_cor <- (tmat_aor$p1+tmat_aor$p2)/2

			#
			try_angles <- rep(NA, length(dptp))
			
			# Number of points moved over the plane
			num_points_over <- rep(NA, length(dptp))
			angle <- NA
			
			for(j in 1:length(cpp)){

				# Use point at max distance to plane to set transformation
				# This will make all point to plane distances either negative or close to negative
				# Faster than trying with all cp points
				cpp_max <- xyz[cpp[j], ]

				# Define circle
				circle <- defineCircle(center=tmat_cor, nvector=tmat_aor$v, point_on_radius=cpp_max)

				# Find intersection of circle and plane - to find vector that will move point into plane about rotational axis
				icp <- intersectCirclePlane(circle=circle, P=xyz[cpv[1],], N=plane_n)
				cpp_max_proj <- circlePoint(circle, icp)

				#if(sum(!is.na(cpp_max_proj[,1])) == 0){
				#	i <- i + 1
				#	next
				#}

				if(sum(!is.na(cpp_max_proj[,1])) == 0) next

				# Select point closest to non-transformed (should be pretty close)
				cpp_max_proj <- cpp_max_proj[which.min(dppt(cpp_max_proj, cpp_max)),]

				# Get angle
				try_angles[j] <- avec(cpp_max-tmat_cor, cpp_max_proj-tmat_cor, axis=tmat_aor$v, about.axis=TRUE)

				# Set transformation matrices
				tmat1 <- tmat2 <- tmat3 <- diag(4)

				# Find rotation transformation
				tmat2[1:3, 1:3] <- tMatrixEP(v=tmat_aor$v, try_angles[j])

				# Get translation transformations
				tmat1[1:3,4] <- tmat_cor
				tmat3[1:3,4] <- -tmat_cor

				# Get transformation matrix
				tmat <- tmat1 %*% tmat2 %*% tmat3

				# Find distance from points to plane (minus means on the side of the normal vector)
				dptp2 <- distPointToPlane(applyTransform(xyz[cpp,], tmat), n=plane_n, q=xyz[cpv[1],])
				
				num_points_over[j] <- sum(dptp2 > -1e-4)
				
				if(num_points_over[j] == length(cpp)){
					angle <- try_angles[j]
					break
				}
				#cat(paste0(j, ') ', try_angles[j], ': ', paste0(dptp2 > 1e-3, collapse=','), '\n'))
				#cat('\n')
			}
			
			if(is.na(angle)){

				if(!any(!is.na(num_points_over))){
					i <- i + 1
					next
				}
				
				angle <- try_angles[which.max(num_points_over)[1]]
			}
			
			#
			dptp_mat[i, 4] <- angle*(180/pi)

			# Set transformation matrices
			tmat1 <- tmat2 <- tmat3 <- diag(4)

			# Find rotation transformation
			tmat2[1:3, 1:3] <- tMatrixEP(v=tmat_aor$v, angle)

			# Get translation transformations
			tmat1[1:3,4] <- tmat_cor
			tmat3[1:3,4] <- -tmat_cor

			# Get transformation matrix
			tmat <- tmat1 %*% tmat2 %*% tmat3

			# Find distance from points to plane (minus means on the side of the normal vector)
			dptp2 <- distPointToPlane(applyTransform(xyz[cpp,], tmat), n=plane_n, q=xyz[cpv[1],])
			dptp_mat[i, 5:7] <- c(min(dptp2, na.rm=TRUE), mean(dptp2, na.rm=TRUE), max(dptp2, na.rm=TRUE))

		}else{
		
			stop("Transformation with constraint plane on dependent (child) body not yet supported.")
		}
		
		i <- i + 1
	}
	
	list(
		'tmat'=tmat,
		'dist'=dptp_mat
	)
}