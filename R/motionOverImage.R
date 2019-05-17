motionOverImage <- function(motion, images, camera, images.save = NULL, landmarks = NULL, 
	meshes = NULL, frame.match = TRUE, scale = 1, animate.speed = 1, body.col = NULL, 
	mesh.opacity = 1, xyz.radius = 1, xyz.opacity = 1, xyz.col = NULL, focal = 150, cone = TRUE, 
	bg.col = 'white', image.opacity = 1, close = TRUE){

	if(!file.exists(images[1])) stop(paste0('Images input file/directory "', images[1],'" not found.'))
	if(!file.exists(meshes[1])) stop(paste0('Meshes input file/directory "', meshes[1],'" not found.'))

	# Get image filenames
	if(file.info(images[1])$isdir){
		image_filenames <- list.files(images[1])
		image_files <- paste0(images[1], '/', image_filenames)
	}else{
		image_filenames <- basename(images)
		image_files <- images
	}

	# Get mesh filenames
	if(!is.null(meshes)){
		if(file.info(meshes[1])$isdir){
			mesh_files <- list.files(meshes[1])
		}else{
			mesh_files <- meshes
		}
	}

	# Read if filepath
	if(length(landmarks) == 1) landmarks <- readLandmarks(landmarks)

	# Match frames between motion and image names using regular expression
	if(frame.match){

		image_frames <- rep(NA, length(image_filenames))
		frames_greg <- gregexpr('[0-9]+', image_filenames)
		for(j in 1:length(image_filenames)){
			image_frames[j] <- as.numeric(substr(image_filenames[j], tail(frames_greg[[j]],1), 
				tail(frames_greg[[j]],1)+tail(attr(frames_greg[[j]],"match.length"),1)-1))
		}

		# Find overlapping frames between read motion and images
		motion_frames_m <- which(motion$frame %in% image_frames)
		image_frames_m <- which(image_frames %in% motion$frame[motion_frames_m])

	}else{

		if(length(motion$frame) != length(image_frames)) stop("If 'frame.match' is TRUE, the number of input images (", length(image_frames), ") must be equal to the number of frames in the input motion object (", length(motion$frame), ").")

		motion_frames_m <- 1:length(motion$frame)
		image_frames_m <- 1:length(image_frames)
	}

	# Check that overlapping frames were found
	if(length(image_frames_m) == 0) warning("No overlapping frames found between motion and images.")

	# Read camera
	camera <- readCam(camera, scale=scale)

	# Set save image names
	if(!is.null(images.save)){

		if(!file.exists(images.save[1])) stop(paste0("images.save input '", images.save[1], "' not found."))

		if(file.info(images.save[1])$isdir){

			# Set save as filenames
			save_images_filenames <- basename(image_files[image_frames_m])

			# Create full filepath
			save_images_as <- paste0(images.save, '/', save_images_filenames)

		}else{

			save_images_as <- images.save

			# Check that lengths match
			if(length(save_images_as) != length(image_frames_m)) stop(paste0("Length of input save as image names (", length(save_images_as), ") does not match number of image frames matching motion frames (", length(image_frames_m), ")."))
		}
	}

	# Get tmat body names
	tmat_bodies <- dimnames(motion$tmat)[[3]]

	# Set body colors
	if(is.null(body.col)) body.col <- setNames(svg.pal(length(tmat_bodies)), tmat_bodies)

	# Make sure number of colors matches
	if(length(body.col) != length(tmat_bodies)){
		if(length(body.col) == 1){
			body.col <- setNames(rep(body.col, length(tmat_bodies)), tmat_bodies)
		}else{
			stop(paste0('The number of body colors (', length(body.col), ') does not match the number bodies in motion$tmat (', 
				length(tmat_bodies), ').'))
		}
	}
	
	# Make sure names match
	tmat_bodies_in <- tmat_bodies %in% names(body.col)
	if(sum(!tmat_bodies_in) > 0) stop(paste0("The following body/bodies was/were not found in names of body.col: '", paste0(tmat_bodies[!tmat_bodies_in], collapse="', '"), "'."))

	# Find matching bodies from mesh filenames
	if(!is.null(meshes)){

		mesh_filenames <- basename(mesh_files)
		meshes_body <- rep(NA, length(mesh_files))
		for(i in 1:length(mesh_filenames)){
			for(j in 1:length(tmat_bodies)){
				if(grepl(paste0('(^| )', tmat_bodies[j], '( |[.])'), mesh_filenames[i])) meshes_body[i] <- tmat_bodies[j]
			}
		}

		# Remove meshes that don't have a matching body in tmat
		mesh_files <- mesh_files[!is.na(meshes_body)]
		meshes_body <- meshes_body[!is.na(meshes_body)]
	}
	
	# Find matching bodies from landmark names
	if(!is.null(landmarks)){
		
		# Get landmark names
		landmark_names <- rownames(landmarks)

		# Get body associated with each landmark
		landmark_body <- unlist(lapply(strsplit(landmark_names, '(_|-)'), head, 1))
		
		# Set rows to keep
		rows_keep <- landmark_body %in% tmat_bodies

		# Limit to landmarks matching body in tmat
		landmarks_m <- landmarks[rows_keep, ]
		landmark_body <- landmark_body[rows_keep]
		
		# Set landmark colors
		if(is.null(xyz.col)){
			xyz.col <- body.col[landmark_body]
		}else{
			xyz.col <- xyz.col[rows_keep]
		}
	}

	# Set distance of plane from pinhole based on max distance from pinhole to points
	xyz_pin_dist <- apply(motion$xyz[,,motion_frames_m], 3, dppt, camera$pinhole)
	cam_dist <- max(xyz_pin_dist, na.rm=TRUE)

	# Start viewer
	#svg.new(file=paste0(study_dir_git, '/Unify animations/', fname_no_ext, '.html'), animate.speed=0.2, mode='webgl')

	if(!is.null(images.save)){
		svg.new(file=save_images_as, mode='webgl', col=bg.col)
	}else{
		svg.new(file=NULL, animate.speed=animate.speed, mode='webgl', col=bg.col)
	}

	# Plot all meshes
	if(!is.null(meshes)){
		for(i in 1:length(mesh_files)){
			svg.mesh(file=mesh_files[i], col=body.col[meshes_body[i]], opacity=mesh.opacity, name=meshes_body[i])
		}
	}

	# Plot landmarks as spheres
	if(!is.null(landmarks)) svg.spheres(x=landmarks_m, name=landmark_body, radius=xyz.radius, col=xyz.col, opacity=xyz.opacity)

	# Transform objects
	svg.transform(tmarr=motion$tmat[,,,motion_frames_m], applyto=dimnames(motion$tmat)[[3]], time=motion$time[motion_frames_m])

	# Add camera
#	svg.camera(camera=camera, focal=focal, image=image_files[image_frames_m], cone=cone, image.opacity=image.opacity, 
#		times=motion$time[motion_frames_m], set=TRUE, plane.dist=cam_dist)

	if(close) svg.close()
}