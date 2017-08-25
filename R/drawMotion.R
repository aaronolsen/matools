drawMotion <- function(file, pts, pts.path, pts.cex = 1, pts.col = 'black', path.col = 'black', 
	path.lwd = 1, frame = TRUE, window.title = 'Draw Motion'){

	svg.new(file, window.title=window.title)

	# Get paths to connect
	read_lines <- readLines(pts.path)
	plist_names <- list()
	for(line_num in 1:length(read_lines)) plist_names[[length(plist_names)+1]] <- strsplit(read_lines[line_num], ',[ ]?')[[1]]

	# Get xyz limits
	ranges <- apply(pts, 2, 'range', na.rm=TRUE)

	svg.points(pts, col.stroke=pts.col, col.fill=pts.col, cex=pts.cex)

	# Draw model paths
	for(i in 1:length(plist_names)){

		path_names <- plist_names[[i]][plist_names[[i]] %in% dimnames(pts)[[1]]]
		if(length(path_names) > 0){

			d <- c()
			for(j in 1:length(path_names)) d <- c(d, which(path_names[j] == dimnames(pts)[[1]]))
			
			svg.pathsC(path=d, col.stroke=path.col, lwd=path.lwd)
		}
	}

	# Draw frame
	if(frame) svg.frame(ranges, z.index=-1)
	
	svg.close()
}