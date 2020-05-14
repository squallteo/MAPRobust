#remember to cite github https://github.com/TGuillerme
bhatt.coeff = function(x,y, bw=bw.nrd0, ...) {
  #SANITIZING
  #x
  if(class(x) != 'numeric') {
    stop("'x' must be numeric.")
  }
  if(length(x) < 2) {
    stop("'x' need at least two data points.")
  }
  
  #y
  if(class(y) != 'numeric') {
    stop("'y' must be numeric.")
  }
  if(length(y) < 2) {
    stop("'y' need at least two data points.")
  }
  
  #bw
  if(length(bw) != 1) {
    stop("'bw' must be either a single numeric value or a single function.")   
  }
  if(class(bw) != 'function') {
    if(class(bw) != 'numeric') {
      stop("'bw' must be either a single numeric value or a single function.")   
    }
  }
  #Avoiding non-entire numbers
  if(class(bw) == 'numeric') {
    bw<-round(bw)
  }
  
  #BHATTACHARYYA COEFFICIENT
  #sum(sqrt(x relative counts in bin_i * y relative counts in bin_i))
  
  #Setting the right number of bins (i)
  if(class(bw) == 'function') {
    #Bin width
    band.width<-bw(c(x,y), ...)
    #Bin breaks
    bin.breaks<-seq(from=min(c(x,y)), to=max(c(x,y)+band.width), by=band.width) #adding an extra bandwith to the max to be sure to include all the data
    #Number of bins
    bin.n<-length(bin.breaks)-1
  } else {
    #Bin breaks
    bin.breaks<-hist(c(x,y), breaks=bw, plot=F)$breaks
    #Bin width
    band.width<-diff(bin.breaks)[1]
    #Number of bins
    bin.n<-bw
  }
  
  #Counting the number of elements per bin
  histx<-hist(x, breaks=bin.breaks, plot=FALSE)[[2]]
  histy<-hist(y, breaks=bin.breaks, plot=FALSE)[[2]]
  #Relative counts
  rel.histx<-histx/sum(histx)
  rel.histy<-histy/sum(histy)
  
  #Calculating the Bhattacharyya Coefficient (sum of the square root of the multiple of the relative counts of both distributions)
  bhatt.coeff<-sum(sqrt(rel.histx*rel.histy))
  return(bhatt.coeff)
  #End
}
