procrustes <- function(x, scale=TRUE, scale.to.mean.Csize=FALSE){

	k <- dim(x)[1]
	m <- dim(x)[2]
	n <- dim(x)[3]
	rownames_x <- rownames(x)

	all_common <- rowSums(!apply(x[, 1, ], 2, is.na)) == n
	#print(names(all_common)[all_common])

	# CENTER
	for(i in 1:n){
		x[,, i] <- x[,, i] - matrix(apply(x[all_common,, i], 2, mean, na.rm=TRUE), k, m, byrow=TRUE)
		#print(apply(x[all_common,, i], 2, mean, na.rm=TRUE)) #Centroid
	}
	
	Csize <- setNames(centroidSize(x[all_common, , ]), dimnames(x)[[3]])
	
	# SCALE
	if(scale){
		for(i in 1:n){
			x[,, i] <- x[,, i] / matrix(sqrt(sum(x[all_common,, i]^2, na.rm=TRUE)), k, m)
			#print(sqrt(sum(x[all_common,, i]^2, na.rm=TRUE))) #Centroid-size
		}
	}
	
	# ROTATE, ALIGN ALL CONFIGURATIONS TO FIRST CONFIGURATION
	for(i in 2:n){
		common <- is.na(x[,1, 1])+is.na(x[,1, i]) == 0
		SVD <- svd(t(x[common,, 1]) %*% x[common,, i])
		L <- diag(SVD$d)
		S <- ifelse(L<0, -1, L)
		S <- ifelse(L>0, 1, L)
		R <- SVD$v %*% S %*% t(SVD$u) # the rotation matrix
		x[,, i] <- x[,, i] %*% R
		SSE <- sum((x[,, i] - x[,, 1])^2, na.rm=TRUE)
		#print(SSE) # Sum of squared differences between target and reference configuration
	}
	mean_config <- t(apply(x, 1, function(x) apply(x, 1, mean)))
	means_SSE <- sum((mean_config - x[,, 1])^2, na.rm=TRUE)

	# ALIGN ALL CONFIGURATIONS TO MEAN CONFIGURATION, RE-ALIGN UNTIL SSE BETWEEN CONSECUTIVE MEAN CONFIGURATIONS < 1e-7
	j <- 0
	while(means_SSE > 1e-7 && j < 5) { # runs a max of 5 iterations or until SSE between means is < 1e-7
		mean_config1 <- t(apply(x, 1, function(x) apply(x, 1, mean)))
		for(i in 1:n){
			common <- is.na(mean_config1[, 1])+is.na(x[,1, i]) == 0
			SVD <- svd(t(mean_config1[common, ]) %*% x[common,, i])
			L <- diag(SVD$d)
			S <- ifelse(L<0, -1, L)
			S <- ifelse(L>0, 1, L)
			R <- SVD$v %*% S %*% t(SVD$u) # the rotation matrix
			x[,, i] <- x[,, i] %*% R
		}
		mean_config2 <- t(apply(x, 1, function(x) apply(x, 1, mean)))
		means_SSE <- sum((mean_config2 - mean_config1)^2, na.rm=TRUE)
		j <- j + 1
	}
	
	if(scale.to.mean.Csize){
		mean_Csize <- mean(Csize, na.rm=TRUE)
		for(i in 1:n) x[,, i] <- x[,, i] * mean_Csize
	}

	return(list('coords'=x, 'Csize'=Csize, 'common'=all_common, 'consensus'=apply(x, 2, 'rowMeans', na.rm=TRUE)))
}