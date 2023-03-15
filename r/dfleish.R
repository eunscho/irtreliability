dfleish <- function(x, mean = 0, sd = 1, skew = 0, kurt = 0) {
  fleishtarget <- function(bcd, skewkurt) {
    b <- bcd[1]
    cc <- bcd[2]
    d <- bcd[3]
    skew <- skewkurt[1]
    kurt <- skewkurt[2]
    out <- (2 - ( 2*b^2 + 12*b*d + skew^2/(b^2+24*b*d+105*d^2+2)^2 + 30*d^2 ) )^2 +
      (kurt - ( 24*(b*d+cc^2*(1+b^2+28*b*d)+d^2*(12+48*b*d+141*cc^2+225*d^2)) ) )^2+
      (cc - (skew/(2*(b^2+24*b*d+105*d^2+2)) ) )^2
    out
  }
  findbcd <- function(skew, kurt){
    optim(c(1, 0, 0), fleishtarget, skewkurt = c(skew,kurt), method = "BFGS",
          control = list(ndeps = rep(1e-10, 3), reltol = 1e-10, maxit = 1e8))
  }
  #####################################################
  #Fleishman Density
  #####################################################
  if(skew == 0 & kurt == 0) {
    out <- dnorm(x, mean = mean, sd = sd)
    } else {
      fcs=findbcd(skew, kurt)$par
      p=-fcs[2]/(3*fcs[3])
      q=p^3+(fcs[2]*fcs[1]-3*fcs[3]*(-fcs[2]-(x-mean)/sd))/(6*fcs[3]^2)
      r=fcs[1]/(3*fcs[3])
      z=(q+sqrt(abs(q^2+(r-p^2)^3)))^(1/3)-abs(q-sqrt(abs(q^2+(r-p^2)^3)))^(1/3)+p
      z0=z
      out <- 1/sd*dnorm(z0)/(fcs[1]+2*fcs[2]*z0+3*fcs[3]*z0^2)
    }
  out
}

