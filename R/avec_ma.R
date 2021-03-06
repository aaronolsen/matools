avec_ma <- function(u, v, axis=NULL, about.axis=FALSE, max.pi=FALSE){

	if(anyNA(u) || anyNA(v)) return(NA)

	if(sqrt(sum(u * u)) == 0) stop("Input vector 'u' is zero-length")
	if(sqrt(sum(v * v)) == 0) stop("Input vector 'v' is zero-length")

	uu <- uvector_ma(u)
	vu <- uvector_ma(v)
	
	c <- sum(uu*vu) / sqrt(sum(uu*uu)) * sqrt(sum(vu*vu))

	if(round(abs(c), digits=12) == 1){
		angle <- 0

		# Check for opposite direction
		if(sum(abs(c(uu[1]+vu[1], uu[2]+vu[2], uu[3]+vu[3]))) < 1e-9){
			angle <- pi
		}

	}else{
		if(max.pi){
			angle <- min(acos(c), pi-acos(c))
		}else{
			angle <- acos(c)
		}
	}

	if(!is.null(axis)){
	
		if(sqrt(sum(axis * axis)) == 0){
			if(abs(angle) < 1e-10) return(0)
			stop("Input vector 'axis' is zero-length")
		}

		axis <- uvector_ma(axis)
		
		if(about.axis){

			# DETERMINE ANGLE BETWEEN VECTORS ABOUT AXIS
			up <- pointPlaneProj_ma(uu, c(0,0,0), axis)
			vp <- pointPlaneProj_ma(vu, c(0,0,0), axis)
			
			if(sqrt(sum(up * up)) == 0 || sqrt(sum(vp * vp)) == 0) return(0)

			return(avec_ma(up, vp, axis=axis, about.axis=FALSE))

		}else{

			# DETERMINE DIRECTION USING AXIS
			if(dppt(uvector_ma(cprod_ma(uu, vu)), axis) < dppt(uvector_ma(cprod_ma(vu, uu)), axis)){
				return(-angle)
			}else{
				return(angle)
			}
		}
	}else{

		# DETERMINE DIRECTION USING CROSS-PRODUCT VECTOR AND EULER 
		if(abs(angle) > 0){
		
			cprod_uv <- cprod_ma(uu, vu)
			um <- sqrt(sum(u^2))
			vm <- sqrt(sum(u^2))
			
			if(dppt((uu %*% tMatrixEP(cprod_uv, angle))*vm, v) <= dppt((uu %*% tMatrixEP(cprod_uv, -angle))*vm, v)){
				return(angle)
			}else{
				return(-angle)
			}
		}
	}
	
	return(angle)
}