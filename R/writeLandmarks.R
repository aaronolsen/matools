writeLandmarks <- function(x, file){

	# Convert array into matrix
	write_mat <- arr2mat(x, dim.to.row=1)
	
	if(grepl('[.]csv$', file)) write.csv(x=write_mat, file=file)
}