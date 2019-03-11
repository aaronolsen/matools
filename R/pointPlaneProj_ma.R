pointPlaneProj_ma <- function(q, p, n){

	if(is.matrix(q)){
		q_ret <- q
		for(i in 1:nrow(q)){
			if(is.na(q[i,1])) next
			q_ret[i,] <- pointPlaneProj_ma(q[i,], p=p, n=n)
		}
		return(q_ret)
	}

	if(q[1] == "non-coincident"){
		warning("Input point non-coincident")
		return(c(NA, NA, NA))
	}
	
	# MAKE SURE N IS A UNIT VECTOR
	n <- uvector_ma(n)

	# PROJECTION OF POINT Q ONTO PLANE DEFINED BY POINT P AND NORMAL VECTOR N
	q - sum((q - p) * n) * n
}