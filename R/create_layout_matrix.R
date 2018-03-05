create_layout_matrix <- function(num.elem, elem.plots, ncol){

	nrow <- ceiling(num.elem / (ncol/elem.plots))

	layout_mat <- matrix(0, nrow, ncol)

	# Fill layout matrix
	rowi <- 1
	coli <- 1
	n <- 1
	for(i in 1:num.elem){
		#cat(rowi, coli, '\n')
		layout_mat[rowi, coli:(coli+elem.plots-1)] <- n:(n+elem.plots-1)
		coli <- coli + elem.plots
		if(coli > ncol){
			coli <- 1
			rowi <- rowi + 1
		}
		n <- n + elem.plots
	}
	
	layout_mat
}
