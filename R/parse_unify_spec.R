parse_unify_spec <- function(ulist, xr_arr, ct_mat, regexp = FALSE, print.progress = FALSE){

	# Check that ulist has element names
	if(is.null(names(ulist))) stop("Unify specification list names are NULL. Input must be a list with named elements in which each name corresponds to an animated element.")

	# Get names of markers from CT mat
	ct_names <- rownames(ct_mat)

	# Create parent-child matrix
	parent_children <- matrix(NA, nrow=0, ncol=2)

	# Make sure ulist is properly formatted
	for(elem_name in names(ulist)){

		# Set nameless single list as 'align'
		if(is.null(names(ulist[[elem_name]])) && length(ulist[[elem_name]]) > 1) stop('Element "', elem_name, '" has multiple sub elements but they are not named (e.g. "align", "point", "orient").')

		# Set nameless single list as 'align'
		if(is.null(names(ulist[[elem_name]]))) ulist[[elem_name]] <- list('align'=ulist[[elem_name]])
		
		for(ltype in c('align', 'point', 'plane')){

			if(is.null(ulist[[elem_name]][[ltype]])) next
			
			if(ltype == 'plane'){
			
				# Check that plane list has two elements
				if(length(ulist[[elem_name]][[ltype]]) != 2) stop('Plane list for "', elem_name, '" has length of ', length(ulist[[elem_name]][[ltype]]), '. Length must be 2.')

				# Check that at least one of objects has same name as element
				if(!elem_name %in% names(ulist[[elem_name]][[ltype]])) stop('One of the two names of the plane list, "', paste0(names(ulist[[elem_name]][[ltype]]), collapse='" and "'), '", must be the same as the corresponding parent element, "', elem_name, '"')
			}
		}
	}

	# Find matches based on regexp matching
	if(regexp){

		for(elem_name in names(ulist)){

			for(l_type in names(ulist[[elem_name]])){
			
				if(l_type == 'align'){

					# Start vector for matching names
					names_match <- c()

					for(i in 1:length(ulist[[elem_name]][[l_type]])){

						# Find match
						which_match <- which(grepl(ulist[[elem_name]][[l_type]][i], ct_names))
					
						if(length(which_match) == 0) warning(paste0('No matching markers found for regular expression: "', ulist[[elem_name]][[l_type]][i], '"'))

						names_match <- c(names_match, ct_names[which_match])
					}
				
					ulist[[elem_name]][[l_type]] <- names_match

				}else{

					for(sub_elem in names(ulist[[elem_name]][[l_type]])){
					
						# Start vector for matching names
						names_match <- c()

						for(i in 1:length(ulist[[elem_name]][[l_type]][[sub_elem]])){
							
							# Find match
							which_match <- which(grepl(ulist[[elem_name]][[l_type]][[sub_elem]][i], ct_names))
							
							if(length(which_match) == 0) warning(paste0('No matching markers found for regular expression: "', ulist[[elem_name]][[l_type]][[sub_elem]][i],'"'))
							
							names_match <- c(names_match, ct_names[which_match])
						}
						
						# Add to list
						ulist[[elem_name]][[l_type]][[sub_elem]] <- names_match

						# Skip if parent
						if(sub_elem == elem_name) next

						# Check whether all points are motion markers
						all_mo_markers <- all(names_match %in% dimnames(xr_arr)[[1]])

						# If markers are all in motion set parent-child doesn't need to be added because markers do not need to be animated with parent
						if(!all_mo_markers) parent_children <- rbind(parent_children, c(sub_elem, elem_name))
					}
				}
			}
		}
	}
	
	for(elem_name in names(ulist)){

		if(length(ulist[[elem_name]][['align']]) == 1){

			if(is.null(ulist[[elem_name]][['point']])){
				stop(paste0('Only a single align marker "', ulist[[elem_name]][['align']], '" and no point-to markers have been assigned to body "', elem_name, '". If only one align marker is assigned then at least 2 point-to markers are needed for unification.\n'))
			}else{
				
				point_markers <- unlist(ulist[[elem_name]][['point']])
				
				if(is.null(ulist[[elem_name]][['plane']]) && length(point_markers) == 1){
					stop(paste0('Only a single align marker "', ulist[[elem_name]][['align']], '" and a single point-to marker "', point_markers, '" have been assigned to body "', elem_name, '". This will result in an underdetermined unification. It would be better to assign at least 2 point-to markers.\n'))
				}
			}
		}

		# Add points to transform with each body (excluding any motion markers)
		if(!is.null(ulist[[elem_name]][['point']])){
			for(sub_elem in names(ulist[[elem_name]][['point']])){
			
				# Get markers to add to transform
				add_to <- ulist[[elem_name]][['point']][[sub_elem]]
				
				# Remove any that are in motion array
				add_to <- add_to[!add_to %in% dimnames(xr_arr)[[1]]]
				
				# Add if non-zero length
				if(length(add_to) > 0) ulist[[sub_elem]][['transform']] <- c(ulist[[sub_elem]][['transform']], add_to)
			}
		}

		if(!is.null(ulist[[elem_name]][['plane']])){
			
			sub_elem <- names(ulist[[elem_name]][['plane']])[names(ulist[[elem_name]][['plane']]) != elem_name]
			
			# Get markers to add to transform
			add_to <- ulist[[elem_name]][['plane']][[sub_elem]]
			
			# Remove any that are in motion array
			add_to <- add_to[!add_to %in% dimnames(xr_arr)[[1]]]
			
			# Add if non-zero length
			if(length(add_to) > 0) ulist[[sub_elem]][['transform']] <- c(ulist[[sub_elem]][['transform']], add_to)
		}
	}
	
	## Set unification order
	if(nrow(parent_children) > 0){

		# Sort by left column - just for readability
		parent_children <- parent_children[order(parent_children[,1]),, drop=FALSE]

		# Set unique body names
		bodies_unique <- names(ulist)
	
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
		
					# Find all parents
					parents_vec <- parent_children[bodies_dep[j] == parent_children[,2], 1]
			
					# Skip if only parent
					if(bodies_dep[j] %in% parents_only) next

					# Skip if all have not been given an order yet
					if(any(is.na(body_dep_order[parents_vec]))) next

					#cat(bodies_dep[j], '\n')
					#cat('\t', paste0(parents_vec, collapse=','), '\n')

					# If order is assigned to parent, assign later order to child
					for(i in 1:length(parents_vec)){
				
						# Get order
						parent_order <- body_dep_order[parents_vec[i]]
				
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
	
		#print(body_dep_order)
	
		# Set order
		body_order <- setNames(c(bodies_indep, names(sort(body_dep_order))), NULL)
	
		# Re-order ulist
		ulist <- ulist[body_order]	
	}
		
	ulist
}
