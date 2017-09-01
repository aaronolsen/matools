addMotion <- function(add, to, add.dim = 1){

	# Check if null
	if(is.null(add)) return(to)
	if(is.null(to)) return(add)

	# Add xyz
	if('xyz' %in% class(add)){
		
		if(is.null(to$xyz)){
		
			# If null, create new
			to$xyz <- add
		
		}else{

			# Set new dimensions
			dim1_names <- unique(c(dimnames(to$xyz)[[1]], dimnames(add)[[1]]))

			if(add.dim == 1){

				# Create new array
				new_xyz <- array(NA, dim=c(length(dim1_names), dim(to$xyz)[2:3]), dimnames=list(dim1_names, 
					dimnames(to$xyz)[[2]], dimnames(to$xyz)[[3]]))
				class(new_xyz) <- 'xyz'

				# Add rows
				new_xyz[dimnames(to$xyz)[[1]], , ] <- to$xyz
				new_xyz[dimnames(add)[[1]], , ] <- add

			}else if(add.dim == 3){
			
				# Create new array
				new_xyz <- array(NA, dim=c(length(dim1_names), dim(to$xyz)[2], dim(to$xyz)[3]+dim(add)[3]), dimnames=list(dim1_names, 
					dimnames(to$xyz)[[2]], c(dimnames(to$xyz)[[3]], dimnames(add)[[3]])))
				class(new_xyz) <- 'xyz'

				# Add rows
				new_xyz[dimnames(to$xyz)[[1]], , 1:dim(to$xyz)[3]] <- to$xyz
				new_xyz[dimnames(add)[[1]], , (dim(to$xyz)[3]+1):(dim(to$xyz)[3]+dim(add)[3])] <- add
			}

			# Add to structure
			to$xyz <- new_xyz
		}

	}else if('tmat' %in% class(add)){

		if(is.null(to$tmat)){
		
			# If null, create new
			to$tmat <- add

		}else{

			# Set new dimensions
			dim3_names <- unique(c(dimnames(to$tmat)[[3]], dimnames(add)[[3]]))

			if(add.dim == 4){
			
				# Create new array
				new_tmat <- array(NA, dim=c(4,4,length(dim3_names),dim(add)[4]+dim(to$tmat)[4]), dimnames=list(NULL, NULL, dim3_names, NULL))
				class(new_tmat) <- 'tmat'

				# If not null, add rows
				new_tmat[, , dimnames(to$tmat)[[3]], 1:dim(to$tmat)[4]] <- to$tmat
				new_tmat[, , dimnames(add)[[3]], (dim(to$tmat)[4]+1):(dim(to$tmat)[4]+dim(add)[4])] <- add
			}

			# Add to structure
			to$tmat <- new_tmat
		}

	}else if('motion' %in% class(add)){

		if(!'motion' %in% class(to)){
			
			# Create new motion structure from add
			if(!is.null(add$xyz)) class(add$xyz) <- 'xyz'
			if(!is.null(add$tmat)) class(add$tmat) <- 'tmat'
			
			return(add)
		}
		
		for(xn in names(add)){
		
			if(xn %in% c('replace.rows', 'n.iter')) next
		
			if('xyz' %in% class(add[[xn]]) || xn == 'xyz'){
				class(add[[xn]]) <- 'xyz'
				to <- addMotion(add[[xn]], to, add.dim=3)
			}else if('tmat' %in% class(add[[xn]]) || xn == 'tmat'){
				class(add[[xn]]) <- 'tmat'
				to <- addMotion(add[[xn]], to, add.dim=4)
			}else{
				to[[xn]] <- c(to[[xn]], add[[xn]])
			}
		}
	}

	to
}
