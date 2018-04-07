bin_replace <- function(x, bin_mat){

	xr <- x
	
	for(rown in 1:nrow(bin_mat)){
		xr[(x >= bin_mat[rown, 1])*(x < bin_mat[rown, 2]) == 1] <- bin_mat[rown, 3]
	}
	
	xr
}