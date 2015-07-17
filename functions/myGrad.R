#################
# FUNCTION INFO #
#################
# PURPOSE: Uses the predicted values of "the.gam" to calculate partial derivatives of x and y (finite differences method)
# AUTHOR: Scott Large 2013
# REVIEWED BY:
# vERSION: 0.1
#

# the.gam <- gam1
my.grad <- function(the.gam, xs, ys, eps.mat) { 
  # Creates an array of predicted values from "the.gam" along x and y sequences "xs", and "ys", for
  # each row in the esp.mat. For partial derivatives of x and y use eps.mat listed below: # 
  # eps <- 1e-7
  # eps.list <- matrix(nrow = 9, ncol = 2, c(+eps, +eps, 0, 0, -eps, 0, -eps, +eps, -eps,
  #                                          +eps, 0, +eps, 0, 0, -eps, -eps, -eps, +eps), byrow = F)
  ### CREATE EPS TERMS ###  
  m.terms <- attr(terms(the.gam), "term.labels") # Select model terms i.e. "X" and "Y" or "Env" and "Dri"
  lp.terms <- the.gam$coefficients # coefficients in the model, to create an array
  p.array <- array(NA, dim = c(sp.len, length(lp.terms), nrow(eps.list))) # sp.len X #coeffecients X eps terms array
  for(i in 1:nrow(eps.mat)) {
    new.dd <- expand.grid(list(x = xs + eps.mat[i,1], # Creates an expanded matrix based on the x and y sequences
                               y = ys + eps.mat[i,2]))
    colnames(new.dd) <- m.terms # Adds the appropriate column names for the gam model
    p.i <- predict.gam(the.gam, newdata = new.dd, type = "lpmatrix") # predicts using the gam for the provided data
    p.array[,,i] <- p.i  # Puts each eps term in the array
  }
  # Create list to hold the data
  nt <- length(m.terms)
  lD <- vector(mode = "list", length = nt + 3)
  names(lD) <- c(m.terms, "xx", "yy", "xy")
  ### CALCULATE Fx, Fx, Fxx, Fyy, and Fxy ###
  Xp <- (p.array[,,2] - p.array[,,5]) / (2*eps) # f(x+h, y) - f(x-h, y)/2h
  Yp <- (p.array[,,3] - p.array[,,6]) / (2*eps) # f(x, y+k) - f(x, y-k)/2k
  XXp <-(p.array[,,2] - 2*(p.array[,,4]) +  p.array[,,5]) / (eps^2) # f(x+h, y) - 2f(x,y) + f(x-h, y) / h^2
  YYp <-(p.array[,,3] - 2*(p.array[,,4]) +  p.array[,,6]) / (eps^2) # f(x, y+h) - 2f(x,y) + f(x, y-h) / k^2
  # f(x+h, y+k) - f(x+h, y) - f(x, y+k) + 2f(x,y) - f(x-h, y) - f(x, y-k) + f(x-h, y-k))/(2hk)
  XYp <- (p.array[,,1] - p.array[,,2] -p.array[,,3] + (2*p.array[,,4]) - p.array[,,5] - p.array[,,6] + p.array[,,7])/ (2*eps*eps) 

  ### MATRIX MULTIPLY THE COEF ###
  Xi <- Xp * 0
  x.want <- grep(m.terms[1], names(lp.terms))
  Xi[, x.want] <- Xp[, x.want]
  x.df <- Xi %*% coef(the.gam)
  x.df.sd <- rowSums(Xi %*% the.gam$Vp * Xi)^.5
  lD[[1]] <- list(grad = x.df, se.grad = x.df.sd)
  # and for Y's
  Yi <- Yp * 0
  y.want <- grep(m.terms[2], names(lp.terms))
  Yi[, y.want] <- Yp[, y.want]
  y.df <- Yi %*% coef(the.gam)
  y.df.sd <- rowSums(Yi %*% the.gam$Vp * Yi)^.5
  lD[[2]] <- list(grad = y.df, se.grad = y.df.sd)
  # # and for Fxx
  XXi <- XXp * 0
  xx.want <- grep(m.terms[1], names(lp.terms))
  XXi[, xx.want] <- XXp[, xx.want]
  xx.df <- XXi %*% coef(the.gam)
  xx.df.sd <- rowSums(XXi %*% the.gam$Vp * XXi)^.5
  lD[[3]] <- list(grad = xx.df, se.grad = xx.df.sd)
  # # and for Fyy
  YYi <- YYp * 0
  yy.want <- grep(m.terms[2], names(lp.terms))
  YYi[, yy.want] <- YYp[, yy.want]
  yy.df <- YYi %*% coef(the.gam)
  yy.df.sd <- rowSums(YYi %*% the.gam$Vp * YYi)^.5
  lD[[4]] <- list(grad = yy.df, se.grad = yy.df.sd)
  # # and for Fxy
  XYi <- XYp * 0
  XYi[, -1] <- XYp[,-1]
  xy.df <- XYi %*% coef(the.gam)
  xy.df.sd <- rowSums(XYi %*% the.gam$Vp * XYi)^.5
  lD[[5]] <- list(grad = xy.df, se.grad = xy.df.sd)
  # Package it all up
  class(lD) <- "my.grad"
  #   lD$gamModel <- the.gam
  #   lD$eps <- eps
  #   lD$eval <- expand.grid(list(x = xs, y = ys))
  return(lD)
}
