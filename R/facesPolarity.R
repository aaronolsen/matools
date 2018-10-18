facesPolarity <- function(t1, t2, return.code = TRUE){

	## Eventually, for each pair of triangles
	#  0 : Triangles intersect
	#  1 : All points on positive side
	# -1 : All points on negative side
	#  2 : Triangles do not intersect but not all points are on one side

	# Get normal vectors
	n1 <- cprod(t1[2,]-t1[1,], t1[3,]-t1[1,])
	n2 <- cprod(t2[2,]-t2[1,], t2[3,]-t2[1,])

	# Get triangle centers
	tc1 <- colMeans(t1)
	tc2 <- colMeans(t2)

	# Check distance from triangle 2's vertices to the triangle 1 plane
	dpp <- rep(NA, 3)
	for(i in 1:3) dpp[i] <- distPointToPlane(t2[i, ], n1, tc1)
	if(sum(dpp < 1e-10) > 0) on_edge <- TRUE

	# Set polarity, default is that t2 crosses t1 plane - but triangles are not intersecting
	polarity <- 2

	# All points are on positive side
	if(sum(dpp > 0) == 3) polarity <- 1

	# All points are on negative side
	if(sum(dpp < 0) == 3) polarity <- -1

	if(polarity != 2){
		if(return.code){ return(polarity) }else{ return(list('polarity'=polarity)) }
	}

	# Check whether triangles intersect
	int_tris <- intersectTriangles(t1, t2, return.logical=return.code)	#, return.logical=TRUE
	#print(int_tris)
	
	if(return.code){
		if(int_tris) return(0)
		return(polarity)
	}

	#
	int_tris
}