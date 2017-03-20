applyTransform <- function(pmat, tarr, assoc){

	# Create array for transformed points
	parr <- array(NA, dim=c(nrow(pmat), 3, dim(tarr)[3]), dimnames=list(rownames(pmat), NULL, NULL))
	
	# Apply transformations for each body
	for(body in 1:dim(tarr)[4]){
		
		# Get body name
		body_name <- dimnames(tarr)[[4]][body]

		# Find points associated with body
		body_assoc <- which(assoc == body_name)
		
		# Skip if no points associated with body
		if(length(body_assoc) == 0) next

		# Get point coordinates as matrix for transformation - coerce to matrix if single point
		pcoor <- rbind(matrix(t(pmat[body_assoc, ]), ncol=length(body_assoc)), rep(1, length(body_assoc)))
		colnames(pcoor) <- rownames(pmat)[body_assoc]
		
	#	pcoor <- pcoor - matrix(c(rowMeans(pcoor[1:3, ]), 0), nrow=nrow(pcoor), ncol=ncol(pcoor))
		
#print(rowMeans(pcoor))
#		break
		
		# Apply transformation
		tcoor <- apply(tarr[, , , body], 3, '%*%', pcoor)

		# Convert to array
		tcoor_arr <- array(tcoor, dim=c(4, length(body_assoc), dim(tarr)[3]))

		# Swap first two dimensions (transpose each "matrix" within array) and remove 1s
		parr[colnames(pcoor), , ] <- aperm(tcoor_arr[1:3, , ], perm=c(2,1,3))
	}

	parr
}
