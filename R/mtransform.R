mtransform <- function(mat, tmat){
	
	if(!is.matrix(mat)) mat <- matrix(mat, nrow=1, ncol=length(mat))

	#pmat <- matrix(1, nrow(mat), ncol(mat)+1)
	#pmat[, 1:3] <- mat
	#mat <- t(tmat %*% t(pmat))[, 1:3]

	# Transform using C++ function
	return_mat <- transform_pts(mat, tmat)

	dimnames(return_mat) <- dimnames(mat)
	
	return_mat
}