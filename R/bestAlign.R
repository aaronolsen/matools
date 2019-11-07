bestAlign <- function(m1, m2, m3 = NULL, sign = NULL, pc.sub=c(0.1, 10000), pt.sub=c(0.5, 500), 
	pt.align = TRUE, print.progress = FALSE, as.cloud = FALSE, reflect = FALSE){

	## See fit_joint_model for revised svd code that fixes bug-- I think I've updated this, based on procAlign

	if(as.cloud){
		m1 <- list('vertices'=m1)
		m2 <- list('vertices'=m2)
	}

	## Mesh alignment
	# Check if m1 is mesh/object
	if(is.list(m1)){

		# Print progress
		if(print.progress) cat(paste0('bestAlign()\n'))
	
		# Check n parameters for consistency
		if(length(pc.sub) == 1){
			if(pc.sub < 0){
				pc.sub <- c(-1, -1)
			}else if(pc.sub > 0 && pc.sub < 1){
				pc.sub <- c(pc.sub, 10000)
			}else{
				pc.sub <- c(0.5, pc.sub)
			}
		}
		if(length(pt.sub) == 1){
			if(pt.sub < 0){
				pt.sub <- c(-1, -1)
			}else if(pt.sub > 0 && pt.sub < 1){
				pt.sub <- c(pt.sub, 500)
			}else{
				pt.sub <- c(0.05, pt.sub)
			}
		}
		
		#	pt.sub: The number of points subsampled from all vertices

		# Reflect, if true
		if(reflect) m2$vertices <- mtransform(m2$vertices, -diag(4))

		# Do an initial alignment using PC axes
		# Get points from vertices		
		m1_pts <- m1$vertices
		m2_pts <- m2$vertices
		
		# Get subset for PCA
		if(pc.sub[1] > -1){
			m1_pc_sub <- m1_pts[seq(1, nrow(m1_pts), length=min(nrow(m1_pts), round(pc.sub[1]*nrow(m1_pts)), pc.sub[2])), ]
			m2_pc_sub <- m2_pts[seq(1, nrow(m2_pts), length=min(nrow(m2_pts), round(pc.sub[1]*nrow(m2_pts)), pc.sub[2])), ]
		}else{
			m1_pc_sub <- m1_pts
			m2_pc_sub <- m2_pts
		}

		# Find centroids		
		m1_centroid <- colMeans(m1_pc_sub)
		m2_centroid <- colMeans(m2_pc_sub)

		# Align by PC axes
		# Perform PCA on vertices
		m1_vectors <- prAxes(m1_pc_sub)$vectors
		m2_vectors <- prAxes(m2_pc_sub)$vectors

		# Set all possible vector orientations
		vec_signs <- list(c(1,1,1), c(-1,1,1), c(1,-1,1), c(1,1,-1), c(-1,-1,1), c(1,-1,-1), c(-1,-1,-1), c(-1,1,-1))

		# Create points from PC for alignment
		m1_pr_pts <- rbind(m1_centroid, m1_vectors + matrix(m1_centroid, 3, 3, byrow=TRUE))
		rownames(m1_pr_pts) <- c('centroid', 'v1', 'v2', 'v3')

		# Try each combination
		pc_errors <- c()
		tmats <- list()
		n <- 1
		for(i in 1:length(vec_signs)){

			# Create points from PC for alignment
			m2_pr_pts <- rbind(m2_centroid, vec_signs[[i]]*m2_vectors + matrix(m2_centroid, 3, 3, byrow=TRUE))
			rownames(m2_pr_pts) <- c('centroid', 'v1', 'v2', 'v3')

			# Find transformation of 2nd to 1st mesh based on PR axes
			tmat_pr <- bestAlign(m1_pr_pts, m2_pr_pts)$tmat
		
			for(angle in c(0, pi/2)){

				# Try rotating about 3rd axis
				# As the most minor axis its orientation is noisier and therefore may not be 
				#  homologously oriented between two similar meshes. This tries rotations of 
				#  90 degrees to see if it produces a lower error
				if(angle != 0){
					tmat_1 <- tmat_2 <- tmat_3 <- diag(4)
					tmat_1[1:3, 4] <- colMeans(m2_pr_pts)
					tmat_2[1:3, 1:3] <- tMatrixEP_ma(vec_signs[[i]][3]*m2_vectors[3,], angle)
					tmat_3[1:3, 4] <- -colMeans(m2_pr_pts)
					tmat_pr_a <- tmat_pr %*% tmat_1 %*% tmat_2 %*% tmat_3
				}else{
					tmat_pr_a <- tmat_pr
				}
		
				# Apply transformation to m2_pc_sub
				m2_pc_sub_pr <- applyTransform(to=m2_pc_sub, tmat=tmat_pr_a)
				
				# Save tmat
				tmats[[n]] <- tmat_pr_a

				# Find error
				pc_errors <- c(pc_errors, min_pt_cloud_dist(m1_pc_sub, m2_pc_sub_pr))
				
				# Advance count
				n <- n + 1
			}
		}
		
		# Save error
		error <- min(pc_errors)

		# Print progress
		if(print.progress) cat(paste0('\tInitial alignment error per point: ', round(error/nrow(m1_pc_sub), 3), '\n'))
		
		# Set final transformation as that which gives lowest error
		tmat_pc <- tmats[[which.min(pc_errors)]]

		# Use point clouds to find best overlap
		if(pt.align){

			# Sub-sample the original vertices
			if(pt.sub[1] > -1){
				m1_sub <- m1_pts[seq(1, nrow(m1_pts), length=min(nrow(m1_pts), round(pt.sub[1]*nrow(m1_pts)), pt.sub[2])), ]
				m2_sub <- m2_pts[seq(1, nrow(m2_pts), length=min(nrow(m2_pts), round(pt.sub[1]*nrow(m2_pts)), pt.sub[2])), ]
			}else{
				m1_sub <- m1_pts
				m2_sub <- m2_pts
			}
			
			# Print number of points used
			if(print.progress) cat(paste0('\tUsing ', nrow(m1_sub), ' points (of ', nrow(m1_pts), ') for point cloud alignment\n'))
			
			# Transform m2 using tmat_pc
			# ***
			m2_sub_pr <- applyTransform(to=m2_sub, tmat=tmat_pc)

			# Set initial parameters
			p_start <- c(0,0,0,0,0,0)

			# Set scale of translation bound
			t_bound <- max(apply(apply(m1_sub, 2, 'range'), 2, 'diff'))*0.2
			
			# Get centroid
			m2_sub_pr_centroid <- colMeans(m2_sub_pr)

			# Find initial error
			#rotate_error_init <- align_pt_cloud_error(p=p_start, m1=m1_sub, m2=m2_sub_pr, center=m2_sub_pr_centroid)
			#print(rotate_error_init)
			
			# Optimize points by rotating 3 axes about real marker
			rotation_fit <- tryCatch(
				expr={
					nlminb(start=p_start, objective=align_pt_cloud_error, m1=m1_sub, m2=m2_sub_pr,
					center=m2_sub_pr_centroid,
					lower=c(p_start[1:3]-t_bound, rep(-2*pi, 3)), upper=c(p_start[1:3]+t_bound, rep(2*pi, 3)))
				},
				error=function(cond) {print(cond);return(NULL)},
				warning=function(cond) {print(cond);return(NULL)}
			)

			# Update error
			error <- rotation_fit$objective

			# Print progress
			if(print.progress) cat(paste0('\tError after point cloud optimization per point: ', round(error/nrow(m1_sub), 3), '\n'))
			
			# Get transformation
			tmat_1 <- tmat_2 <- tmat_3 <- diag(4)
			tmat_1[1:3, 4] <- m2_sub_pr_centroid+rotation_fit$par[1:3]
			tmat_2[1:3, 1:3] <- rotationMatrixZYX_ma(rotation_fit$par[4:6])
			tmat_3[1:3, 4] <- -m2_sub_pr_centroid
			tmat <- tmat_1 %*% tmat_2 %*% tmat_3 %*% tmat_pc

		}else{
			tmat <- tmat_pc
		}

		# If m2 points were reflected add reflection before other transformations
		if(reflect){
			reflect_tmat <- -diag(4)
			reflect_tmat[4,4] <- 1
			tmat <- tmat %*% reflect_tmat
		}

		rlist <- list(
			'error'=error,
			'tmat'=tmat
		)
		return(rlist)
	}

	# If m1 is 3-d array
	if(length(dim(m1)) == 3 && length(dim(m2)) == 2){

		# Get number of iterations
		n_iter <- dim(m1)[3]

		# Get common points
		common_names <- rownames(m1)[rownames(m1) %in% dimnames(m2)[[1]]]
		
		# Array for transformed m2
		m2r <- array(NA, dim=c(dim(m2), n_iter), dimnames=list(dimnames(m2)[[1]], NULL, NULL))
		
		# Create transformation array
		tmat <- array(NA, dim=c(4,4,n_iter))

		# Not finished
		stop('Not finished writing')
	}

	if(length(dim(m2)) == 3){
		
		# Get number of iterations
		n_iter <- dim(m2)[3]
		
		# Get common points
		common_names <- rownames(m1)[rownames(m1) %in% dimnames(m2)[[1]]]

		#
		rlist <- list(
			mat=array(NA, dim=dim(m2), dimnames=dimnames(m2)),
			pos.errors=array(NA, dim=c(dim(m1[common_names, ]), n_iter), dimnames=list(common_names, NULL, NULL)),
			dist.errors=array(NA, dim=c(dim(m1[common_names, ])[1], 1, n_iter), dimnames=list(common_names, NULL, NULL)),
			tmat=array(NA, dim=c(4,4,n_iter))
		)

		# Each iteration
		for(iter in 1:n_iter){

			# Find translation and rotation to align m2 to m1
			#m2[, , iter] <- bestAlign(m1[common_names, ], m2[common_names, , iter], m2[, , iter])$mc
			best_align <- bestAlign(m1[common_names, ], m2[common_names, , iter], m2[, , iter])
			
			# Save results
			rlist$mat[,,iter] <- best_align$mc
			rlist$pos.errors[,,iter] <- best_align$pos.errors
			rlist$dist.errors[,,iter] <- best_align$dist.errors
			rlist$tmat[,,iter] <- best_align$tmat
		}

		return(rlist)
	}
	
	# IF INPUTTING A MATRIX WITH ALL ZEROS IN ONE DIMENSION, MAKE IT M2, NOT M1

	if(is.null(m3)) m3r <- NULL

	# SET INITIAL COMMON POINT MATRIX VALUES
	m1o <- m1
	m2o <- m2

	# USE ROWNAMES, IF GIVEN, TO REMOVE NON-CORRESPONDING POINTS
	if(!is.null(rownames(m1)) && !is.null(rownames(m2))){

		# Check for duplicate rownames
		if(length(unique(rownames(m1))) < nrow(m1)) stop('Input parameter "m1" contains duplicate row names.')
		if(length(unique(rownames(m2))) < nrow(m2)) stop('Input parameter "m2" contains duplicate row names.')
		
		# If empty row names, assume same row names as other matrix
		if(any(rownames(m1o) == "")){
			if(!any(rownames(m2o) == "") && nrow(m1o) == nrow(m2o)){
				rownames(m1o) <- rownames(m2o)
			}else{
				stop('Input parameter "m1" contains empty row names and number of rows do not match between m1o and m2o.')
			}
		}
		if(any(rownames(m2o) == "")){
			if(nrow(m1o) == nrow(m2o)){
				rownames(m2o) <- rownames(m1o)
			}else{
				stop('Input parameter "m2" contains empty row names and number of rows do not match between m1o and m2o.')
			}
		}
		
		m1o <- m1o[sort(rownames(m1o)), ]
		m2o <- m2o[sort(rownames(m2o)), ]

		# REMOVE NA VALUES
		m1o <- m1o[!is.na(m1o[, 1]), ]
		m2o <- m2o[!is.na(m2o[, 1]), ]

		m1o[rownames(m1o)[!rownames(m1o) %in% rownames(m2o)], ] <- NA
		m2o[rownames(m2o)[!rownames(m2o) %in% rownames(m1o)], ] <- NA

	}else{

		# REPLACE NON-COMMON LANDMARKS BETWEEN TWO MATRICES WITH NA
		m1o[which(is.na(m2o))] <- NA
		m2o[which(is.na(m1o))] <- NA
	}

	# REMOVE NA VALUES
	m1o <- m1o[!is.na(m1o[, 1]), ]
	m2o <- m2o[!is.na(m2o[, 1]), ]

	# CREATE FIRST TRANSLATION TRANSFORMATION MATRIX
	tmat1 <- diag(4);tmat1[1:3, 4] <- -colMeans(m2o, na.rm=TRUE)
	
	# CENTER M2 ABOUT CENTROID OF COMMON POINTS
	m2c <- mtransform(m2, tmat1)

	# APPLY TRANSLATION TO EXTRA MATRIX, IF PROVIDED
	if(!is.null(m3)) m3c <- mtransform(m3, tmat1)
	
	# CENTER COMMON POINTS
	m1oc <- scale(m1o, center=TRUE, scale=FALSE)
	m2oc <- scale(m2o, center=TRUE, scale=FALSE)
	
	# Find best alignment with just two points
	if(nrow(m1oc) == 2){
		
		# Check if points are identical
		tmat2 <- diag(4)
		if(sum(colMeans(abs(m2oc - m1oc))) > 1e-10){

			m_axis <- cprod_ma(m1oc[2,]-m1oc[1,], m2oc[2,]-m2oc[1,])

			# Find angle between points
			m_avec <- avec_ma(m2oc[2,], m1oc[2,], axis=m_axis, about.axis=TRUE)

			# Find rotation
			RM <- tMatrixEP_ma(v=m_axis, -m_avec)
			
			#rotated <- rotateBody(m2oc, m_axis, -m_avec)
			#cat('-----------\n')
			#print(m1oc)
			#print(rotated)

			# Set rotation matrix			
			tmat2[1:3, 1:3] <- t(RM)
		}

		m2r <- mtransform(m2c, tmat2)

	}else{

		# FIND ROTATION MATRIX TO APPLY TO M2 THAT MINIMIZES DISTANCE BETWEEN M1 AND M2
		SVD <- svd(t(na.omit(m1oc)) %*% na.omit(m2oc))

		# CORRECTION TO ENSURE A RIGHT-HANDED COORDINATE SYSTEM
		S <- diag(3)
		S[3,3] <- sign(det(SVD$v %*% t(SVD$u)))

		# GET ROTATION MATRIX
		# MIGHT CHANGE POINTS RELATIVE TO ONE ANOTHER (VERY SLIGHTLY)
		# I THINK IT ONLY HAPPENS WHEN ONE DIMENSION OF M1 IS ALL ZEROS
		# CAUSES PROBLEM IN DETERMINING FIT PERHAPS
		RM <- SVD$v %*% S %*% t(SVD$u)

		# CREATE ROTATION TRANSFORMATION MATRIX
		tmat2 <- diag(4);tmat2[1:3, 1:3] <- t(RM)

		# TEST ALIGNMENT
		#t2 <- m2c %*% RM
		#print(m1oc[!is.na(m1oc[, 1]), ])
		#print(t2[!is.na(t2[, 1]), ])

		# ROTATE ALL CENTER LANDMARKS IN M2
		m2r <- mtransform(m2c, tmat2)
		if(!is.null(m3)) m3r <- mtransform(m3c, tmat2)

		# TEST WHETHER CHIRALITY OF POINT SET HAS FLIPPED
		if(nrow(m1o) == 3 && sum(!is.na(m2[,1])) > 3){
		
			# IF THE ALIGNMENT FLIPPED THE SET TO ITS MIRROR IMAGE, THE CROSS PRODUCT OF THE 
			#	FIRST THREE POINTS WILL MAINTAIN THE SAME ORIENTATION. BUT ANY OTHER POINTS WILL
			#	BE FLIPPED RELATIVE TO THIS VECTOR. SO IF THE DISTANCE BETWEEN THE CROSS PRODUCT
			#	AND THESE POINTS CHANGES, IT INDICATES THE CHIRALITY HAS BEEN FLIPPED

			# FIND NORMAL VECTORS FOR PRE AND POST ROTATED SETS
			m2c_cprod <- uvector_ma(cprod_ma(m2c[2, ]-m2c[1, ], m2c[3, ]-m2c[1, ]))
			m2r_cprod <- uvector_ma(cprod_ma(m2r[2, ]-m2r[1, ], m2r[3, ]-m2r[1, ]))
		
			# FIND DISTANCE FROM CPROD VECTOR TO OTHER POINTS
			dpp <- dppt(m2c_cprod, m2c[4:min(7,nrow(m2c)), ])
			dpp_r <- dppt(m2r_cprod, m2r[4:min(7,nrow(m2r)), ])
		
			# CHIRALITY HAS FLIPPED, FLIP 3RD COLUMN OF SVD$v AND RE-TRANSFORM
			if(sum(round(abs(dpp - dpp_r), 7)) > 0.001){
				SVD$v[, 3] <- -SVD$v[, 3]
				RM <- SVD$v %*% S %*% t(SVD$u)
				tmat2[1:3, 1:3] <- t(RM)
				m2r <- mtransform(m2c, tmat2)
				if(!is.null(m3)) m3r <- mtransform(m3c, tmat2)
			}else{
			}
		}
	}

	m2or <- mtransform(m2oc, tmat2)

	# CREATE SECOND TRANSLATION TRANSFORMATION MATRIX
	tmat3 <- diag(4);tmat3[1:3, 4] <- colMeans(m1o, na.rm=TRUE)

	# APPLY TRANSLATION PARAMETERS
	m2r <- mtransform(m2r, tmat3)
#	m2r <- m2r + matrix(colMeans(m1o, na.rm=TRUE), nrow=nrow(m2r), ncol=ncol(m2r), byrow=TRUE)
	if(!is.null(m3)) m3r <- mtransform(m3r, tmat3)

	# GET ALIGNMENT ERROR
	errors <- m1oc - m2or
	attr(errors, "scaled:center") <- NULL
	
	dist.errors <- matrix(sqrt(rowSums(errors^2)), ncol=1, dimnames=list(rownames(errors), NULL))

	# GET FINAL TRANSFORMATION MATRIX
	tmat <- tmat3 %*% tmat2 %*% tmat1

	#errors <- m1o - mtransform(m2o, tmat)
	
	list(
		mat=m2r,
		pos.errors=errors,
		dist.errors=dist.errors,
		mc=m3r,
		tmat=tmat
	)
}