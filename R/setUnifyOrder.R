setUnifyOrder <- function(markers){
	
	# Get bodies
	bodies <- unlist(lapply(strsplit(markers, '_|-'), head, 1))

	# Get unique bodies
	bodies_unique <- unique(bodies)

	# Get virtual points
	markers_vp <- markers[grepl('-', markers)]
	
	# Get constraint plane points
	markers_cpv1 <- markers[grepl('_cp[0-9]+v1', markers)]
	markers_cpp1 <- markers[grepl('_cp[0-9]+p1', markers)]
	
	# Get cp IDs
	if(length(markers_cpv1) > 0){
		cpv_ids <- rep(NA, length(markers_cpv1))
		cpp_ids <- rep(NA, length(markers_cpp1))
		for(i in 1:length(markers_cpv1)){
			cpv_ids[i] <- as.numeric(gsub('cp|v[0-9]+[f]?', '', strsplit(markers_cpv1[i], '_')[[1]][2]))
			cpp_ids[i] <- as.numeric(gsub('cp|p[0-9]+[f]?', '', strsplit(markers_cpp1[i], '_')[[1]][2]))
		}
	}else{
		cpv_ids <- NULL
		cp_list <- NULL
	}

	# If no virtual points or constraint planes, any order works; return alphabetic order
	if(length(markers_vp) == 0 && length(markers_cpv1) == 0) return(sort(bodies_unique))

	# Create parent-child matrix
	parent_children <- matrix(NA, nrow=length(markers_vp)+length(markers_cpv1), ncol=2)
	
	# Fill matrix
	# Virtual points
	if(length(markers_vp) > 0){
		for(i in 1:length(markers_vp)) parent_children[i, ] <- strsplit(markers_vp[i], '_|-')[[1]][1:2]
		row_i <- length(markers_vp)
	}else{
		row_i <- 1
	}

	# Constraint plane/points
	if(!is.null(cpv_ids)){
	
		parent_type <- rep(NA, length(cpv_ids))
		child_type <- rep(NA, length(cpv_ids))
		for(i in 1:length(cpv_ids)){
			
			# Find corresponding cp point
			which_cpp <- which(cpp_ids == cpv_ids[i])

			# Find whether vertices or points are fixed (parent)
			parent_type[i] <- ifelse(grepl('f$',markers_cpv1[i]), 'v', 'p')
			child_type[i] <- ifelse(parent_type[i] == 'v', 'p', 'v')

			# Add parents and children
			if(parent_type[i] == 'v'){
				parent_children[i+row_i, ] <- c(strsplit(markers_cpv1[i], '_')[[1]][1], strsplit(markers_cpp1[which_cpp], '_')[[1]][1])
			}else{
				parent_children[i+row_i, ] <- c(strsplit(markers_cpp1[which_cpp], '_')[[1]][1], strsplit(markers_cpv1[i], '_')[[1]][1])
			}
		}

		# Create data frame
		cp_list <- list('parent'=parent_children[(1+row_i):nrow(parent_children),1], 
			'parent.type'=parent_type, 'child'=parent_children[(1+row_i):nrow(parent_children),2], 
			'child.type'=child_type, 'id'=cpv_ids)
		
		#if(length(markers_vp) > 0){
		#	cp_children <- parent_children[(1+row_i):nrow(parent_children),2]
		#	cp_children_in_vp <- cp_children[cp_children %in% parent_children[1:row_i,2]]
		#	if(length(cp_children_in_vp) > 0){
		#		stop(paste0("Body(ies) \"", paste0(cp_children_in_vp, collapse='", "'), "\" are children of both constraint planes and virtual points. Trouble creating hierarchical order of unification."))
		#	}
		#}
	}
	
	# Remove duplicate rows
	rownames(parent_children) <- paste0(parent_children[,1], '-', parent_children[,2])
	parent_children <- parent_children[unique(rownames(parent_children)), ]

	# Get list of all bodies in parent-child matrix
	parent_children_bodies <- unique(c(parent_children))
	
	# Find all independent bodies
	bodies_indep <- bodies_unique[!bodies_unique %in% parent_children_bodies]

	# Find all bodies that are only parents
	parents_only <- unique(parent_children[!parent_children[,1] %in% parent_children[,2], 1])

	## Get children iteratively through parents
	# Set parents only as next parents to start
	next_parents <- parents_only

	n_max <- 100
	n <- 1
	bodies_dep <- c()
	while(n < n_max){
	
		# Get children of parents
		next_parents <- parent_children[parent_children[, 1] %in% next_parents, 2]
		
		# If no children, stop
		if(length(next_parents) == 0) break
		
		# Add to 
		bodies_dep <- c(bodies_dep, next_parents)
		n <- n + 1
	}
	
	list(
		'order'=setNames(c(bodies_indep, parents_only, bodies_dep), NULL),
		'cp'=cp_list
	)
}