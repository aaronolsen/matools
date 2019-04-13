prAxes <- function(pts){

	# PC analysis
	pca <- prcomp(pts)
	
	# Find centroid
	meanX <- apply(pts, 2, mean) 

	# Get points at ends of each axis
	axes <- list()
	ends <- matrix(NA, nrow=2*ncol(pca$x), ncol=length(meanX))
	vectors <- matrix(NA, nrow=ncol(pca$x), ncol=length(meanX))

	j <- 1
	for(i in 1:ncol(pca$x)){
		axes[[i]] <- rbind(meanX + min(pca$x[, i])*pca$rotation[, i], meanX + max(pca$x[, i])*pca$rotation[, i])
		ends[j:(j+1), ] <- axes[[i]]
		vectors[i,] <- uvector_ma(axes[[i]][2,]-axes[[i]][1,])
		j <- j + 2
	}
	
	list('axes'=axes, 'ends'=ends, 'vectors'=vectors)
}