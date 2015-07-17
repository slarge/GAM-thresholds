#################
# FUNCTION INFO #
#################
# PURPOSE: Create a convex hull around raw data with a buffer of 1/diameter 
# AUTHOR: Scott Large 2013
# REVIEWED BY:
# vERSION: 0.1
#

hullFun <- function(hull.df, diameter = 1, npoints = 10) {
  
  circleFun <- function(center = c(0,0), diameter = 1, npoints = 10){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }
  
  # Scale and center original values
  scaled.df <- data.frame(x = scale(hull.df[,1], scale = T, center = T),
                          y = scale(hull.df[,2], scale = T, center = T))
  
  # Take convex hull of the scaled and centered values
  hpts <- chull(scaled.df)
  hpts <- c(hpts, hpts[1])
  # Create a circle of points surrounding each convex hull point
  circ <- data.frame()
  for(cf in 1:length(hpts)) {
    cf.df <- c(scaled.df[hpts[cf],1], scaled.df[hpts[cf],2])
    circ <- rbind(circ, circleFun(center = cf.df, diameter = diameter, npoints = npoints))
  }
  
  hpts2 <- chull(circ)
  hpts2 <- c(hpts2,hpts2[1])
  hull <- data.frame(x = circ[hpts2,1] * sd(hull.df[,1]) + mean(hull.df[,1]),
                     y = circ[hpts2,2] * sd(hull.df[,2]) + mean(hull.df[,2]))
  return(hull)
}
