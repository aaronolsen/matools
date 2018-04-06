slide_fun <- function(data, window, step, fun = 'mean', fill.na = TRUE){

  # From: http://stats.stackexchange.com/questions/3051/mean-of-a-sliding-window-in-r

  total <- length(data)
  spots <- seq(from=round(window/2), to=(total-round(window/2)), by=step)
  result <- rep(NA, length(data))
  for(i in 1:length(spots)){
    result[spots[i]] <- do.call(fun, list(data[(spots[i]-(window/2)):(spots[i]+(window/2))]))
  }
  
  result[1:(round(window/2)-1)] <- result[!is.na(result)][1]
  result[(length(result)-round(window/2)+1):length(result)] <- result[!is.na(result)][sum(!is.na(result))]
  
  result
}