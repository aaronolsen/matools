sampleSphereSurface <- function(n, radius = 1, center = c(0,0,0)){

	# https://stackoverflow.com/questions/5408276/sampling-uniformly-distributed-random-points-inside-a-spherical-volume

	sample_length <- 3000

	seq2pi <- seq(0,2*pi,length=sample_length)
	seq11 <- seq(-1,1,length=sample_length)
	seq01 <- seq(0,1,length=sample_length)

	# Start with slightly more values because some will be 0-length and become NA
	n_over <- round(n*1.01)+10

	phi2 <- sample(seq2pi, size=n_over, replace=TRUE)
	costheta2 <- sample(seq11, size=n_over, replace=TRUE)
	u2 <- sample(seq01, size=n_over, replace=TRUE)
	theta2 <- acos(costheta2)
	v2 <- u2^(1/3)

	pts2 <- matrix(NA, nrow=n_over, ncol=3)
	pts2[, 1] <- v2 * sin(theta2) * cos(phi2)
	pts2[, 2] <- v2 * sin(theta2) * sin(phi2)
	pts2[, 3] <- v2 * cos(theta2)
	pts2 <- matrix(center, nrow=nrow(pts2), ncol=ncol(pts2), byrow=TRUE) + radius*(pts2 / sqrt(apply(pts2^2, 1, sum)))

	# Remove NAs
	pts2 <- pts2[!is.na(pts2[, 1]), ]

	# Trim to right length
	pts2 <- pts2[1:n, ]

	pts2
}