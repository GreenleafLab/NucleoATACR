
#' mean_smooth
#' Smooth using a flat window
#'
#' @param x vector of values
#' @param window size of window for smoothing
#' @return returns vector of size x that has been smoothed.  For values at edge, mean is for window included in x. 
#' (i.e. for fist value, smoothed valued is mean of that value and the window/2 values after)
#' @export
mean_smooth<-function(x, window = 11){
  w = floor(window/2)
  sapply(1:length(x), function(y) mean(x[max(1,y-w):min(length(x),y+w)], na.rm=T))
}



#' tabulate2
#' Modified version of tabulate so that it gives counts in a given span and can handle negative integers
#'
#' @param x vector of values
#' @param min_val minimum value
#' @param max_val maximum value
#' @return returns vector of length max_val - min_val + 1 with values representing counts for values between min_val and max_val
#' @export
tabulate2<-function(x,min_val,max_val){
  if (max_val <= min_val){
    stop("max_val must be greater than min_val")
  }
  if (min_val<0 && max_val >0){
    n = rev(tabulate(-1*(x))[1:(-min_val)])
    p = tabulate(x)[1:max_val]
    z = length(which(x == 0))
    out = c(n,z,p)
    names(out)=min_val:max_val
    return(out)}
  else if (min_val==0 && max_val >0){
    p = tabulate(x)[1:max_val]
    z = length(which(x == 0))
    out = c(z,p)
    names(out)=min_val:max_val
    return(out)}
  else if (min_val > 0 && max_val >0){
    p = tabulate(x)[min_val:max_val]
  }
  else if (min_val <0 && max_val == 0){
    n = rev(tabulate(-1*(x))[1:(-min_val)])
    z = length(which(x == 0))
    out = c(n,z)
    names(out)=min_val:max_val
    return(out)}
  else if (min_val <0 && max_val < 0){
    n = rev(tabulate(-1*(x))[1:(-min_val)])
    out = n
    names(out)=min_val:max_val
    return(out)}
  else{
    stop("something may be amiss with min_val or max_val")
  }
} 