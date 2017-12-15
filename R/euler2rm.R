euler2rm <- function(x, y = NULL, z = NULL, order = 'xyz'){

	if(length(x) == 1){
		x[2] <- y
		x[3] <- z
	}

	z <- matrix(c(cos(x[1]), -sin(x[1]), 0, sin(x[1]), cos(x[1]), 0, 0, 0, 1), 3, 3, byrow=TRUE)
	y <- matrix(c(cos(x[2]), 0, sin(x[2]), 0, 1, 0, -sin(x[2]), 0, cos(x[2])), 3, 3, byrow=TRUE)
	x <- matrix(c(1, 0, 0, 0, cos(x[3]), -sin(x[3]), 0, sin(x[3]), cos(x[3])), 3, 3, byrow=TRUE)

	if(order == 'xyz') return(z %*% y %*% x)
	if(order == 'zyx') return(x %*% y %*% z)
}