sampleOVec <- function(v, n = NULL){

	if(is.null(n)) if(is.matrix(v)){ n <- nrow(v) }else{ n <- 1 }

	# If v is single vector, expand to same length as n_over
	if(is.matrix(v)){
		if(nrow(v) == 1){
			v <- matrix(c(v), nrow=n, 3, byrow=TRUE)
		}else if(ncol(v) == 1){
			v <- matrix(c(v), nrow=n, 3, byrow=TRUE)
		}else{
			if(nrow(v) != n) stop('Number of rows in v (', nrow(v), ') is not equal to n (', n, ').')
		}
	}

	if(!is.matrix(v)){
		vec_input <- TRUE
		v <- matrix(v, nrow=n, 3, byrow=TRUE)
	}else{
		vec_input <- FALSE
	}

	# Get sample of unit vectors
	if(n == 1){
		u_vecs <- matrix(sampleSphereSurface(n=n), nrow=1)
	}else{
		u_vecs <- sampleSphereSurface(n=n)
		#u_vecs[2, ] <- v[2, ] # Test parallel vector detection
	}

	# Create orthogonal vector matrix
	o_vecs <- matrix(NA, nrow(u_vecs), 3)

	# Get cross product
	for(i in 1:nrow(o_vecs)) o_vecs[i, ] <- uvector(cprod(v[i,], u_vecs[i,]))

	niter <- 1
	while(sum(rowSums(abs(o_vecs)) == 0) > 0 || niter < 10){

		# Remove zero rows (vector pairs that were parallel)
		for(i in 1:nrow(o_vecs)){
			if(sum(o_vecs[i, ]) != 0) next
			o_vecs[i, ] <- sampleOVec(1, v[i,])
		}

		niter <- niter + 1
	}

	if(vec_input && nrow(o_vecs) == 1) o_vecs <- c(o_vecs)
	o_vecs
}