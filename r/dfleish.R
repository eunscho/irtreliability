fleishtarget<-function(x,a){
  b<-x[1];cc<-x[2];d<-x[3];g1<-a[1];g2<-a[2]
  (2 - ( 2*b^2 + 12*b*d + g1^2/(b^2+24*b*d+105*d^2+2)^2 + 30*d^2 ) )^2 +
    (g2 - ( 24*(b*d+cc^2*(1+b^2+28*b*d)+d^2*(12+48*b*d+141*cc^2+225*d^2)) ) )^2+
    (cc - (g1/(2*(b^2+24*b*d+105*d^2+2)) ) )^2
}
findbcd<-function(skew,kurtosis){
  optim(c(1,0,0),fleishtarget,a=c(skew,kurtosis),method="BFGS",
        control=list(ndeps=rep(1e-10,3),reltol=1e-10,maxit=1e8))
}
fcs=findbcd(skew,kurt)$par
#####################################################
#Fleishman Density
#####################################################
dfleish=function(x,skew,kurt,mt,st){
  if(skew==0 & kurt==0) {dnorm(x,mean=mt,sd=st)} else
    if(skew!=0 | kurt!=0){
      fcs=findbcd(skew,kurt)$par
      p=-fcs[2]/(3*fcs[3])
      q=p^3+(fcs[2]*fcs[1]-3*fcs[3]*(-fcs[2]-(x-mt)/st))/(6*fcs[3]^2)
      r=fcs[1]/(3*fcs[3])
      z=(q+sqrt(abs(q^2+(r-p^2)^3)))^(1/3)-abs(q-sqrt(abs(q^2+(r-p^2)^3)))^(1/3)+p
      z0=z
      1/st*dnorm(z0)/(fcs[1]+2*fcs[2]*z0+3*fcs[3]*z0^2)
    }}
