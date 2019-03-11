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
	if(length(markers_vp) == 0 && length(markers_cpv1) == 0) return(list('order'=sort(bodies_unique)))

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
	if(!is.matrix(parent_children)) parent_children <- matrix(parent_children, 1, 2)

	# Get list of all bodies in parent-child matrix
	parent_children_bodies <- unique(c(parent_children))
	
	# Find all independent bodies
	bodies_indep <- bodies_unique[!bodies_unique %in% parent_children_bodies]

	# Find all bodies that are only parents
	parents_only <- unique(parent_children[!parent_children[,1] %in% parent_children[,2], 1])

	# Find all dependent bodies
	bodies_dep <- bodies_unique[bodies_unique %in% parent_children_bodies]

	# Set order of dependent bodies
	body_dep_order <- setNames(rep(NA, length(bodies_dep)), bodies_dep)

	# Set parents as zero
	body_dep_order[parents_only] <- 0

	#print(parent_children)
	#print(body_dep_order)
	
	# Set order of dependent bodies
	# Loop until all are assigned
	n <- 0
	while(n < 10){

		if(length(body_dep_order) > 1){
			for(j in 1:length(body_dep_order)){
		
				#cat(bodies_dep[j], '\n')

				# Find all parents
				parents_vec <- parent_children[bodies_dep[j] == parent_children[,2], 1]
			
				#cat('\t', paste0(parents_vec, collapse=','), '\n')

				# Skip if only parent
				if(bodies_dep[j] %in% parents_only) next
			
				# If order is assigned to parent, assign later order to child
				for(i in 1:length(parents_vec)){
				
					# Get order
					parent_order <- body_dep_order[parents_vec[i]]
				
					# Skip if not given an order yet (NA)
					if(is.na(parent_order)) next
				
					#cat('\t', parent_order, '\n')
				
					if(is.na(body_dep_order[bodies_dep[j]])){
						body_dep_order[bodies_dep[j]] <- parent_order+1
					}else{
						if(parent_order+1 > body_dep_order[bodies_dep[j]]) body_dep_order[bodies_dep[j]] <- parent_order+1
					}
				}

				#cat('-------\n')
			}
		}
		
		if(!any(is.na(body_dep_order))) break
		
		n <- n + 1
	}
		
	list(
		'order'=setNames(c(bodies_indep, names(sort(body_dep_order))), NULL),
		'cp'=cp_list
	)
}