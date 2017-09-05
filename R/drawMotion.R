drawMotion <- function(file, pts, pts.path = NULL, pts.cex = 1, pts.col.stroke = 'black', 
	pts.col.fill = 'black', path.col = 'black', path.lwd = 1, pts2 = NULL, pts2.cex = 1, 
	pts2.col.stroke = 'red', pts2.col.fill = 'red', frame = TRUE, duration = 1, draw.labels = FALSE, 
	window.title = 'Draw Motion'){

	svg.new(file, animate.duration=duration, window.title=window.title)

	# Get xyz limits
	ranges <- apply(pts, 2, 'range', na.rm=TRUE)

	# Draw first set of points
	svg.points(pts, col.stroke=pts.col.stroke, col.fill=pts.col.fill, cex=pts.cex)

	# Draw model paths
	if(!is.null(pts.path)){

		read_lines <- readLines(pts.path)
		plist_names <- list()
		for(line_num in 1:length(read_lines)) plist_names[[length(plist_names)+1]] <- strsplit(read_lines[line_num], ',[ ]?')[[1]]
	
		for(i in 1:length(plist_names)){

			path_names <- plist_names[[i]][plist_names[[i]] %in% dimnames(pts)[[1]]]
			if(length(path_names) > 0){

				d <- c()
				for(j in 1:length(path_names)) d <- c(d, which(path_names[j] == dimnames(pts)[[1]]))
			
				svg.pathsC(path=d, col.stroke=path.col, lwd=path.lwd)
			}
		}
	}

	# Draw labels
	if(draw.labels) if(!is.null(dimnames(pts)[[1]])) svg.text(pts, labels=dimnames(pts)[[1]], font.size=1, col='blue')

	# Draw second set of points
	if(!is.null(pts2)) svg.points(pts2, col.stroke=pts2.col.stroke, col.fill=pts2.col.fill, cex=pts2.cex)

	# Draw frame
	if(frame) svg.frame(ranges, z.index=-1)
	
	svg.close()
}