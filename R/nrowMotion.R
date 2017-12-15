nrowMotion <- function(motion){

	# Get number of rows
	if(!is.null(motion$n.iter)){

		num_rows <- motion$n.iter

	}else{

		for(xn in names(motion)){

			if(xn %in% c('replace.rows', 'n.iter')) next

			if(is.null(dim(motion[[xn]]))){
				num_rows <- length(motion[[xn]])
			}else{
				if(length(dim(motion[[xn]])) == 2){
					num_rows <- dim(motion[[xn]])[2]
				}else if(length(dim(motion[[xn]])) == 3){
					num_rows <- dim(motion[[xn]])[3]
				}else if(length(dim(motion[[xn]])) == 4){
					num_rows <- dim(motion[[xn]])[4]
				}
			}
		}
	}

	num_rows
}