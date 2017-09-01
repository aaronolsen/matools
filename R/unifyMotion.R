unifyMotion <- function(ct_mat, xr_arr, add.xr = FALSE, print.progress = TRUE, print.progress.iter = c(1)){

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

	# Set body associations
	body_assoc <- body_names

	# Get unique body names
	body_names <- unique(body_names)
	
	# Create transformation matrix names
	tm_names <- body_names

	# Create transformation matrix array
	tm_arr <- array(NA, dim=c(4, 4, length(tm_names), dim(xr_arr)[3]), 
		dimnames=list(NULL, NULL, tm_names, NULL))

	# Create coordinate system names
	cs_names <- body_names

	# Create array for initial coordinate system positions and orientations
	cs_ini <- array(NA, dim=c(4, dim(xr_arr)[2], length(cs_names)), dimnames=list(NULL, NULL, cs_names))

	# Fill coordinate system matrix - initial configuration
	for(cs_name in dimnames(cs_ini)[[3]]){
	
		if(cs_name %in% body_names){

			# Find all markers associated with bodys
			ct_mat_sub <- ct_mat[grepl(paste0(cs_name, '[_|-]'), paste0(rownames(ct_mat), '_')), ]
			
			# Find body centroid
			if(is.matrix(ct_mat_sub)){
				cs_ini[1, , cs_name] <- colMeans(ct_mat_sub, na.rm=TRUE)
			}else{
				cs_ini[1, , cs_name] <- ct_mat_sub
			}
			
			# Add orientation points
			cs_ini[2:4, , cs_name] <- matrix(cs_ini[1, , cs_name], nrow=3, ncol=3, byrow=TRUE) + diag(3)
		}
	}

	# Create array for transformed CT markers
	ct_arr <- array(NA, dim=c(dim(ct_mat)[1], dim(xr_arr)[2:3]), dimnames=list(rownames(ct_mat), dimnames(xr_arr)[[3]]))

	# Order in which to align body points - do bodies with virtual markers first
	body_order <- setNames(rep(2, length(body_names)), body_names)
	body_order[unique(names(body_names_vm[!is.na(body_names_vm)]))] <- 1
	body_order <- sort(body_order)

	# Get names of virtual markers
	virtual_markers <- rownames(ct_mat)[!is.na(body_names_vm)]

	# Add virtual markers to xr_arr
	xr_arr_n <- array(NA, dim=c(dim(xr_arr)[1]+length(virtual_markers), dim(xr_arr)[2], dim(xr_arr)[3]), 
		dimnames=list(c(dimnames(xr_arr)[[1]], virtual_markers), dimnames(xr_arr)[[2]], dimnames(xr_arr)[[3]]))
	xr_arr_n[dimnames(xr_arr)[[1]], , ] <- xr_arr

	# Unification errors
	errors <- matrix(NA, nrow=dim(xr_arr_n)[3], ncol=length(body_order), dimnames=list(NULL, names(body_order)))

	# Iterate through each frame of video
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

		# Do unification of all common points to approximately orient bodies with single X-ray marker
		m1 <- xr_arr[rownames(ct_mat)[rownames(ct_mat) %in% dimnames(xr_arr)[[1]]], , iter]
		m2 <- ct_mat
		m3 <- apply(cs_ini, 2, as.matrix)

		# Skip if all markers are NA
		if(sum(is.na(m1)) == length(m1)) next

		align_ct_xr <- bestAlign(m1, m2, m3=m3, sign=1)
		ct_mat_align <- align_ct_xr$mat
		cs_ini_align <- cs_ini
		for(rowi in seq(1, nrow(m3), by=4)) cs_ini_align[, , ((rowi-1) / 4) + 1] <- align_ct_xr$mc[rowi:(rowi+3), ]
		
		# Unify markers
		for(body_name in names(body_order)){

			#if(!body_name %in% c('SuspensoriumL', 'Neurocranium')) next
			if(print.progress && iter %in% print.progress.iter) cat(paste0('\t', body_name))
			
			# Get body landmarks
			#ct_mat_sub <- ct_mat[which(body_assoc == body_name), ]
			ct_row_sub <- grepl(paste0(body_name, '[_|-]'), paste0(rownames(ct_mat), '_'))
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
			
			# CT marker names in Xray
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

				# Find translation from CT coordinates to X-ray coordinates
				#tmat <- diag(4)

				tmat <- align_ct_xr$tmat
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
				align <- bestAlign(xr_mat_sub, ct_mat_sub, m3=cs_ini[, , body_name], sign=1)
				ct_mat_sub_t <- align$mat
				
				# Save error
				errors[iter, body_name] <- mean(align$dist.errors)

				# Save transformation matrix
				tm_arr[, , body_name, iter] <- align$tmat
			}

			if(print.progress && iter %in% print.progress.iter && !is.null(align)){
				#print(align$dist.error)
				#cat(paste0('\t\tError range: ', paste(round(range(align$dist.error), 3), collapse=', '), '\n'))
				#cat(paste0('\t\t', paste0(rownames(align$dist.error), ': ', c(round(align$dist.error, 3)), collapse='\n\t\t'), '\n'))
				cat(paste0('\t\tErrors: ', paste0(rownames(align$dist.error), ': ', c(round(align$dist.error, 3)), collapse='; '), '\n'))
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

			ct_arr[rownames(ct_mat_sub_t), , iter] <- ct_mat_sub_t

			#print(xr_mat_sub);print(ct_mat_sub_t)
		}
		
		#return(1)
		
		#if(iter > 4) break
		
		#print(ct_arr[, , iter])
	}

	# If any values in CT array are NA, replace with X-ray markers that have the same name
	ct_is_na <- is.na(ct_arr[, 1, 1])
	if(sum(ct_is_na) > 0){
		ct_na_in_xr <- names(ct_is_na)[ct_is_na][names(ct_is_na)[ct_is_na] %in% dimnames(xr_arr)[[1]]]
		if(length(ct_na_in_xr) > 0) ct_arr[ct_na_in_xr, , ] <- xr_arr[ct_na_in_xr, , ]
	}
	
	# Add X-ray points to animated CT points	
	if(add.xr){
		ct_arr_new <- array(NA, dim=c(dim(ct_arr)[1] + dim(xr_arr_n)[1], dim(ct_arr)[2], dim(ct_arr)[3]),
			dimnames=list(c(dimnames(ct_arr)[[1]],dimnames(xr_arr_n)[[1]]), NULL, NULL))
		ct_arr_new[dimnames(ct_arr)[[1]], , ] <- ct_arr
		ct_arr_new[dimnames(xr_arr_n)[[1]], , ] <- xr_arr_n
		ct_arr <- ct_arr_new
	}

	list(
		'xyz.ct'=ct_arr,
		'xyz.xr'=xr_arr_n,
		'tmat'=tm_arr,
		'error'=errors
	)
}