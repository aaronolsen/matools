unifyMotion <- function(motion, xyz.mat, print.progress = TRUE, print.progress.iter = c(1), 
	replace.xyz = TRUE, plot.diag = NULL, vp.use = TRUE, cp.use = TRUE, skip.bodies = c(),
	near.bodies = NULL){

	add.xr <- TRUE

	# Set point array
	xr_arr <- motion$xyz
	
	# Set matrix coordinates
	ct_mat <- xyz.mat

	# Remove CT coordinates that are NA
	ct_mat <- ct_mat[!is.na(ct_mat[, 1]), ]

	# Get body names
	body_names <- rownames(ct_mat)
	body_names <- gsub('_[A-Za-z0-9-]*', '', body_names)

	#body_names <- gsub('[0-9]', '', body_names)
	#for(i in 1:2) body_names <- gsub(paste0('_(ant|sup|mid|inf|pos)[_]?'), '_', body_names)
	#body_names <- gsub('_$', '', body_names)

	# Get virtual markers
	body_names_vm <- rep(NA, length(body_names))
	body_names_vm[grepl('-', body_names)] <- gsub('[A-Za-z]+-', '', body_names[grepl('-', body_names)])
	body_names <- gsub('-[A-Za-z]+', '', body_names)
	names(body_names_vm) <- body_names

	# If vp.use is FALSE, remove virtual markers
	if(!vp.use){
		ct_mat <- ct_mat[is.na(body_names_vm),]
		body_names <- rownames(ct_mat)
		body_names <- gsub('_[A-Za-z0-9-]*', '', body_names)
		body_names_vm <- rep(NA, length(body_names))
		names(body_names_vm) <- body_names
	}
	
	# Remove any skip body names
	if(length(skip.bodies) > 0) body_names <- body_names[!body_names %in% skip.bodies]

	# Set body associations
	body_assoc <- body_names

	# Get unique body names
	body_names <- unique(body_names)
	
	# Create transformation matrix names
	tm_names <- body_names

	# Create transformation matrix array
	tm_arr <- array(NA, dim=c(4, 4, length(tm_names), dim(xr_arr)[3]), 
		dimnames=list(NULL, NULL, tm_names, NULL))

	# Create array for transformed CT markers
	ct_arr <- array(NA, dim=c(dim(ct_mat)[1], dim(xr_arr)[2:3]), dimnames=list(rownames(ct_mat), dimnames(xr_arr)[[3]]))

	# Order in which to align body points - do bodies with virtual markers first
	#body_order <- setNames(rep(2, length(body_names)), body_names)
	#body_order[unique(names(body_names_vm[!is.na(body_names_vm)]))] <- 1
	#body_order <- sort(body_order)
	
	# Set order in which to unify bodies
	set_unify_order <- setUnifyOrder(rownames(ct_mat))
	unify_order <- set_unify_order$order
	unify_cp <- set_unify_order$cp
	
	# If non NULL create array to store results of constraint plane transformations
	if(cp.use && !is.null(unify_cp)){
		col_names <- c('pre.min', 'pre.mean', 'pre.max', 'angle', 'post.min', 'post.mean', 'post.max')
		cp_array <- array(NA, dim=c(length(unify_cp$id), length(col_names), motion$n.iter), dimnames=list(NULL, col_names, NULL))
	}


	# Get names of virtual markers
	virtual_markers <- rownames(ct_mat)[!is.na(body_names_vm)]

	# Add virtual markers to xr_arr
	xr_arr_n <- array(NA, dim=c(dim(xr_arr)[1]+length(virtual_markers), dim(xr_arr)[2], dim(xr_arr)[3]), 
		dimnames=list(c(dimnames(xr_arr)[[1]], virtual_markers), dimnames(xr_arr)[[2]], dimnames(xr_arr)[[3]]))
	xr_arr_n[dimnames(xr_arr)[[1]], , ] <- xr_arr

	# Convert list elements into regular expression match
	if(!is.null(near.bodies)){
		for(i in 1:length(near.bodies)){
			near.bodies[[i]] <- grepl(paste0('^', paste0(near.bodies[[i]], collapse='|'), '(_|-)'), dimnames(xr_arr)[[1]])
		}
	#	near.bodies <- setNames(unlist(near.bodies), names(near.bodies))
	}

	# Unification errors
	errors <- matrix(NA, nrow=dim(xr_arr_n)[3], ncol=length(unify_order), dimnames=list(NULL, unify_order))

	# Iterate through each frame
	min_iter <- 1
	max_iter <- dim(xr_arr_n)[3]
	for(iter in min_iter:max_iter){
	
		if(print.progress){
			if(iter %in% print.progress.iter){
				cat(paste0('Iteration: ', iter, '\n'))
			}else if(iter == max(print.progress.iter) + 1){
				#cat(paste0('Iterations: ', iter))
			}else{
				#cat(iter)
			}
		}

		# Unify markers
		for(body_name in unify_order){

			#if(!body_name %in% c('SuspensoriumL', 'Neurocranium')) next
			if(print.progress && iter %in% print.progress.iter) cat(paste0('\t', body_name))
			
			# Get all landmarks associated with body, including virtual point
			ct_row_sub <- which(grepl(paste0('^', body_name, '(_|-)'), rownames(ct_mat)))
			ct_row_sub <- c(ct_row_sub, which(grepl(paste0('-', body_name), rownames(ct_mat))))
			ct_mat_sub <- ct_mat[ct_row_sub, ]
			#print(rownames(ct_mat)[grepl(paste0(body_name, '_|-'), rownames(ct_mat))])
		
			# Find which CT markers are in XR array
			ct_in_xr <- rownames(ct_mat_sub) %in% dimnames(xr_arr_n)[[1]]

			if(sum(ct_in_xr) == 0){
				if(print.progress && iter %in% print.progress.iter) cat(paste0(': 0 common point(s) between CT and X-Ray sets for body \'', body_name, '\'\n'))
				next
			}

			# If fewer than 3 non-NA values, skip
			#if(sum(!is.na(xr_arr_n[rownames(ct_mat_sub)[ct_in_xr], 1, iter])) < 3){
			#	if(print.progress && iter %in% print.progress.iter) print(rownames(ct_mat_sub)[ct_in_xr])
			#	next
			#}
			
			# CT marker names in Xray (beads and virtual points, does not include constraint planes/points)
			ct_in_xr_names <- rownames(ct_mat_sub)[ct_in_xr]

			# XROMM markers in CT
			xr_mat_sub <- matrix(xr_arr_n[ct_in_xr_names, , iter], nrow=sum(ct_in_xr), ncol=ncol(xr_arr_n), 
				dimnames=list(ct_in_xr_names, NULL))
			
			# Remove NA values
			xr_mat_sub <- matrix(xr_mat_sub[!is.na(xr_mat_sub[, 1]), ], nrow=sum(!is.na(xr_mat_sub[, 1])), ncol=ncol(xr_mat_sub), 
				dimnames=list(rownames(xr_mat_sub)[!is.na(xr_mat_sub[, 1])], NULL))

			# Translate bodies with single marker in common
			if(nrow(xr_mat_sub) == 0){

				if(print.progress && iter %in% print.progress.iter) cat(': all X-ray markers are NA\n')
				
				next

			}else if(nrow(xr_mat_sub) == 1){
			
				# Check that marker is not virtual
				common_marker <- rownames(xr_mat_sub)
				if(grepl('-', common_marker)){
					if(print.progress && iter %in% print.progress.iter) cat(': single marker is virtual marker\n')
					next
				}

				if(print.progress && iter %in% print.progress.iter) cat(': translate body with CT alignment\n')

				if(print.progress && iter %in% print.progress.iter) cat('\t\tSetting body orientation based on alignment with specified neighboring landmarks\n')

				# Do unification of set common points to approximately orient bodies with single X-ray marker
				if(is.null(near.bodies) && is.null(near.bodies[[body_name]])){
					m1 <- xr_arr[rownames(ct_mat)[rownames(ct_mat) %in% dimnames(xr_arr)[[1]]], , iter]
				}else{
					m1 <- xr_arr[near.bodies[[body_name]], , iter]
				}

				m2 <- ct_mat

				# Skip if all markers are NA
				if(sum(is.na(m1)) == length(m1)) next
				
				align_ct_xr <- bestAlign(m1, m2, sign=1) #m3=m3, 
				ct_mat_align <- align_ct_xr$mat

				# Set initial transformation based on alignment
				tmat <- align_ct_xr$tmat

				# Translate based on single marker
				tmat[1:3, 4] <- tmat[1:3, 4] + (xr_mat_sub[common_marker, ] - ct_mat_align[common_marker, ])

				# Translate CT markers based on common point with X-ray markers
				ct_mat_sub_t <- ct_mat_align[rownames(ct_mat_sub), ] + matrix(xr_mat_sub[common_marker, ] - ct_mat_align[common_marker, ], nrow=nrow(ct_mat_sub), ncol=ncol(ct_mat_sub), byrow=TRUE)

				# Save transformation matrix
				tm_arr[, , body_name, iter] <- tmat

				# Set align as NULL
				align <- NULL

			}else if(nrow(xr_mat_sub) == 2){
			
				# Find which are NA
				which_is_na <- is.na(xr_arr_n[ct_in_xr_names, 1, 1])

				if(print.progress && iter %in% print.progress.iter) cat(paste0(': ', sum(ct_in_xr), ' common point(s) between CT and X-Ray sets for body \'', body_name, '\'\n'))
				if(print.progress && iter %in% print.progress.iter && sum(which_is_na) > 0) cat(paste0('\t\t"', paste0(ct_in_xr_names[which_is_na], collapse='","'), '" is/are are NA\n'))
				next

			}else{

				if(print.progress && iter %in% print.progress.iter) cat(': transform CT markers to align with XROMM markers\n')

				# Transform CT markers to correspond with XROMM markers
				align <- bestAlign(xr_mat_sub, ct_mat_sub, sign=1)	#, m3=cs_ini[, , body_name]
				ct_mat_sub_t <- align$mat
				
				# Save error
				errors[iter, body_name] <- mean(align$dist.errors)

				if(body_name == 'HyoidL'){
					#print(xr_mat_sub)
					#errors[iter, body_name] <- align$dist.errors['Basihyal-HyoidL_Basihyal_hypohyal_jt',]
					#errors[iter, body_name] <- align$dist.errors['HyoidL_bead_cau',]
					#errors[iter, body_name] <- align$dist.errors['HyoidL_bead_cra',]
				}

				# Save transformation matrix
				tm_arr[, , body_name, iter] <- align$tmat
			}

			if(print.progress && iter %in% print.progress.iter && !is.null(align)){
				#print(align$dist.error)
				#cat(paste0('\t\tError range: ', paste(round(range(align$dist.error), 3), collapse=', '), '\n'))
				#cat(paste0('\t\t', paste0(rownames(align$dist.error), ': ', c(round(align$dist.error, 3)), collapse='\n\t\t'), '\n'))
				cat(paste0('\t\tErrors: ', paste0(rownames(align$dist.error), ': ', c(round(align$dist.error, 3)), collapse='; '), '\n'))
			}

			#if(print.progress && iter %in% print.progress.iter) print(ct_mat_sub_t)
			
			ct_arr[rownames(ct_mat_sub_t), , iter] <- ct_mat_sub_t

			## Impose any constraint planes
			# Check if body is child in a constraint plane set
			if(cp.use && body_name %in% unify_cp$child){

				# Get CP IDs
				cp_ids <- unify_cp$id[body_name == unify_cp$child]
				
				# Get matching index
				cp_match <- unify_cp$id == cp_ids[1]

				if(print.progress && iter %in% print.progress.iter){
					if(length(cp_ids) == 1){
						cat(paste0('\t\tRefining motion using constraint plane "', cp_ids,'" and parent body "', unify_cp$parent[cp_match], '"\n'))
					}else{
						cat(paste0('\t\tRefining motion using constraint planes "', paste0(cp_ids, collapse=','), '" and parent body "', unify_cp$parent[cp_match], '"\n'))
					}
				}
				
				# Find transformation from constraint plane and points
				cpt <- suppressWarnings(constraintPlaneTransform(ct_arr[, , iter], cp_ids, unify_cp$child.type[cp_match]))
				tmat <- cpt$tmat
				cp_array[cp_ids,,iter] <- cpt$dist
				#print(cp_array[,,iter])

				if(print.progress && iter %in% print.progress.iter) cat(paste0('\t\tApplying transformation to: ', paste0(rownames(ct_mat_sub_t), collapse=', '), '\n'))

				# Apply transformation to points and transformation
				ct_arr[rownames(ct_mat_sub_t), , iter] <- applyTransform(ct_arr[rownames(ct_mat_sub_t), , iter], tmat)
				ct_mat_sub_t <- applyTransform(ct_mat_sub_t, tmat)
				tm_arr[, , body_name, iter] <- tmat %*% align$tmat

				# Get new errors
				beads_and_vp <- rownames(ct_mat_sub_t)[grepl('_bead|-', rownames(ct_mat_sub_t))]
				dist_errors <- distPointToPoint(xr_arr_n[beads_and_vp, , iter], ct_mat_sub_t[beads_and_vp, ])
				if(print.progress && iter %in% print.progress.iter) cat(paste0('\t\tErrors: ', paste0(names(dist_errors), ': ', c(round(dist_errors, 3)), collapse='; '), '\n'))

				# Save new errors
				errors[iter, body_name] <- mean(dist_errors, na.rm=TRUE)
			}

			# Add new positions of any virtual markers to xr_arr_n
			if(sum(rownames(ct_mat_sub) %in% virtual_markers) > 0){

				# Virtual markers in subset
				virtual_markers_sub <- rownames(ct_mat_sub)[rownames(ct_mat_sub) %in% virtual_markers]

				# Remove virtual markers that are non-NA in xr_arr (previously positioned with first body)
				virtual_markers_sub <- virtual_markers_sub[is.na(xr_arr_n[virtual_markers_sub, 1, iter])]

				# Add new virtual markers
				if(length(virtual_markers_sub) > 0){

					if(print.progress && iter %in% print.progress.iter) cat(paste0('\t\tAdding virtual markers: ', paste0(virtual_markers_sub, collapse=', '), '\n'))

					xr_arr_n[virtual_markers_sub, , iter] <- ct_mat_sub_t[virtual_markers_sub, ]
				}
			}

			#print(xr_mat_sub);print(ct_mat_sub_t)
		}
		
		#return(1)
		
		#break
		#if(iter > 4) break
		
		#print(ct_arr[, , iter])
	}

	if(cp.use && !is.null(unify_cp) && print.progress){
		cat(paste0('\nResults of constraint plane transformations:\n'))

		for(i in 1:dim(cp_array)[1]){
		
			cat(paste0('\t', i, ') ', unify_cp$parent[i]), '(')

			if(unify_cp$parent.type[i] == 'v'){ cat('plane') }else{ cat('points') }
			cat(paste0(') - ', unify_cp$child[i], ' ('))
			if(unify_cp$child.type[i] == 'v'){ cat('plane') }else{ cat('points') }
			cat(paste0(')\n'))

			# Skip if all values are NA
			if(sum(!is.na(cp_array[i,,])) == 0){
				cat(paste0('\t\tAll values NA\n'))
				next
			}

			cat(paste0('\t\tPoint-to-plane distance from: Min:', round(min(cp_array[i,'pre.min',], na.rm=TRUE), 3), ', Mean:', round(mean(cp_array[i,'pre.mean',], na.rm=TRUE), 3), ', Max:', round(max(cp_array[i,'pre.max',], na.rm=TRUE), 3), '\n'))
			if(sum(!is.na(cp_array[i,'post.min',])) > 0){
				cat(paste0('\t\tTo: Min:', round(min(cp_array[i,'post.min',], na.rm=TRUE), 3), ', Mean:', round(mean(cp_array[i,'post.mean',], na.rm=TRUE), 3), ', Max:', round(max(cp_array[i,'post.max',], na.rm=TRUE), 3), '\n'))
				cat(paste0('\t\tAngle (deg): Min:', round(min(cp_array[i,'angle',], na.rm=TRUE), 3), ', Mean:', round(mean(cp_array[i,'angle',], na.rm=TRUE), 3), ', Max:', round(max(cp_array[i,'angle',], na.rm=TRUE), 3), '\n'))
			}
			cat(paste0('\t\tNumber of frames corrected: ', sum(!is.na(cp_array[i,'post.min',])), ' of ', motion$n.iter, '\n'))
		}
	}

	# If any values in CT array are NA, replace with X-ray markers that have the same name
	ct_is_na <- is.na(ct_arr[, 1, 1])
	if(sum(ct_is_na) > 0){
		ct_na_in_xr <- names(ct_is_na)[ct_is_na][names(ct_is_na)[ct_is_na] %in% dimnames(xr_arr)[[1]]]
		if(length(ct_na_in_xr) > 0) ct_arr[ct_na_in_xr, , ] <- xr_arr[ct_na_in_xr, , ]
	}
	
	# Add X-ray points to animated CT points	
	if(add.xr){
	
		# Get all marker names
		unique_names <- unique(c(dimnames(ct_arr)[[1]],dimnames(xr_arr_n)[[1]]))
		
		#
		ct_arr_new <- array(NA, dim=c(length(unique_names), dim(ct_arr)[2], dim(ct_arr)[3]),
			dimnames=list(unique_names, NULL, NULL))

		# 
		if(replace.xyz){
			ct_arr_new[dimnames(xr_arr_n)[[1]], , ] <- xr_arr_n
			ct_arr_new[dimnames(ct_arr)[[1]], , ] <- ct_arr
		}else{

			xr_arr_n_non_na <- dimnames(xr_arr_n)[[1]][rowSums(is.na(xr_arr_n[, 1, ])) == 0]

			ct_arr_new[dimnames(ct_arr)[[1]], , ] <- ct_arr
			ct_arr_new[xr_arr_n_non_na, , ] <- xr_arr_n[xr_arr_n_non_na, , ]
		}
		ct_arr <- ct_arr_new
	}

	#
	if(!is.null(plot.diag)){

		# If plot filename doesn't end in pdf, add
		if(!grepl('[.]pdf$', plot.diag, ignore.case=TRUE)) plot.diag <- paste0(plot.diag, '.pdf')

		# Set which columns to plot - not columns that are all NA
		cols_plot <- which(colSums(is.na(errors)) < nrow(errors))

		## Create error diagnostic plot
		pdf(plot.diag, height=length(cols_plot)*2.5, width=max(nrow(errors)/90, 7))
		layout(cbind(1:length(cols_plot)))
		par(mar=c(4.5,4.5,3,1))
	
		for(i in cols_plot){

			if(!any(!is.na(errors[, i]))){
				plot(c(0,1), c(0,1), type='n', main=colnames(errors)[i], 
					xlab='Frame', ylab='Unification error (mm)')
				next
			}

			plot(c(1, nrow(errors)), c(0, max(errors[, i], na.rm=TRUE)), type='n', 
				main=paste0(colnames(errors)[i], ' (Mean: ', round(mean(errors[, i], na.rm=TRUE), 3), '; SD: ', round(sd(errors[, i], na.rm=TRUE), 3), ')'), 
				xlab='Frame', ylab='Unification error (mm)')

			abline(h=0, lty=2, col=gray(0.5))

			points(1:nrow(errors), errors[, i], type='l')
		}
	
		dev.off()
	
	}

	motion[['xyz']] <- ct_arr
	motion[['tmat']] <- tm_arr

	rlist <- list(
		'motion'=motion,
		'error'=errors
		#''=xr_arr_n
	)
}