quat2AxisAngle <- function(q){

	## Converts quaternion to axis of rotation and angle
	
	# Make sure q is normalized
	q <- uvector_ma(q)

	# Get angle
	angle <- 2 * acos(q[4])

	# Get scalar
	scalar <- sqrt(1-q[4]^2)

	# Test to avoid divide by zero, scalar is always positive due to sqrt
	# If scalar close to zero then direction of axis not important
	# If it is important that axis is normalised then replace with x=1 y=z=0
	aor <- c(q[1:3])

	# Scale
	if(scalar > 1e-4) aor <- aor / scalar
	
	list('axis'=aor, 'angle'=angle)
}