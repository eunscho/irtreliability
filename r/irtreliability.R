IRTReliabilityNN=function(J,k,bs,as,c,skew,kurt,mt,st){
  #####################################################
  #Fleishman Functions for Coefficients
  #####################################################
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
  #####################################################
  #Modified Graded Response Model and Derivatives
  #####################################################
  Pijstars=function(x,m,k,a,b,c){
    if(m==0) g=1
    if(m>0 & m<k+1) g=1/(1+exp(-a*(x-(b-c))))
    if(m==k) g=0
    g}
  cs=c(0,c,0)#This allows movement through thresholds in for loop

  dPijstars=function(x,m,k,a,b,c){
    if(m==0) g=0
    if(m>0 & m<k+1) g=a*Pijstars(x,m,k,a,b,c)*(1-Pijstars(x,m,k,a,b,c))
    if(m==k) g=0
    g}

  Pij=function(x,m,k,a,b,c0,c1) Pijstars(x,m-1,k,a,b,c=c0)-Pijstars(x,m,k,a,b,c=c1)
  Aij=function(x,m,k,a,b,c0,c1) Pijstars(x,m-1,k,a,b,c=c0)+Pijstars(x,m,k,a,b,c=c1)
  dAij = function(x,m,k,a,b,c0,c1) dPijstars(x,m-1,k,a,b,c=c0)+dPijstars(x,m,k,a,b,c=c1)
  dPij=function(x,m,k,a,b,c0,c1) a*Pij(x,m,k,a,b,c0,c1)*(1-Aij(x,m,k,a,b,c0,c1))
  d2Pij=function(x,m,k,a,b,c0,c1) a*(dPij(x,m,k,a,b,c0,c1)*(1-Aij(x,m,k,a,b,c0,c1) ) - Pij(x,m,k,a,b,c0,c1)*dAij(x,m,k,a,b,c0,c1))

  #####################################################
  #Information Functions and Reliability of theta hat
  #####################################################
  IIk=function(x,k,a,b,m,c0,c1){
    sapply(x,function(x){
      (dPij(x,m,k,a,b,c0=c0,c1=c1))^2/Pij(x,m,k,a,b,c0=c0,c1=c1) - d2Pij(x,m,k,a,b,c0=c0,c1=c1)
    })
  }

  vthetahatf=function(x,k,as,bs,cs,J,st){
    val=numeric(length(x))
    for(j in 1:J){
      for(i in 1:k){
        val=IIk(x,k=k,a=as[j],b=bs[j],m=i,c0=cs[i],c1=cs[i+1])+val
      }}
    #st^2/(st^2+1/val)
    1/val
  }

  IRTxx=function(x,k,as,bs,cs,J,skew,kurt,mt,st){
    ifelse(x>mt+5*st | x< mt-5*st,0,vthetahatf(x,k,as,bs,cs,J=J,st)*dfleish(x,skew,kurt,mt,st=st))
  }
  vtheta.x=integrate(IRTxx,-Inf,Inf,k,as,bs,cs,J,skew,kurt,mt=mt,st=st)$value
  thetaxx=st^2/(st^2+vtheta.x)
  #####################################################################
  #CTT Reliability
  #####################################################################
  #####################################################################
  #CTT Item Conditional Variance
  #####################################################################
  Eu.theta=function(x,k,a,b,cs){
    sapply(x,function(x) {
      val=numeric(length(x))
      for(m1 in 1:k){
        val=m1*Pij(x=x,m=m1,k,a=a,b=b,c0=cs[m1],c1=cs[m1+1])+val
      }
      val
    })
  }

  CTTIIF=function(x,k,a,b,cs){
    sapply(x,function(x) {
      sval=numeric(length(x))
      for(m1 in 1:k){
        sval=(m1-Eu.theta(x=x,k=k,a=a,b=b,cs=cs))^2*Pij(x=x,m=m1,k,a=a,b=b,c0=cs[m1],c1=cs[m1+1])+sval
      }
      sval
    })
  }
  #####################################################################
  #CTT Test Conditional Standard Error of Measurement
  #####################################################################
  CTTTIF=function(x,k,J,as,bs,cs){
    sapply(x,function(x) {
      fval=0
      for(i in 1:J){
        fval= CTTIIF(x,k,a=as[i],b=bs[i],cs=cs)+fval
      }
      sqrt(fval)
    })
  }


  ECTTTIF=function(x,k,J,as,bs,cs,skew,kurt,mt,st){
    CTTTIF(x,k,J,as,bs,cs)*dfleish(x,skew,kurt,mt,st)
  }

  Vx.t=integrate(ECTTTIF,-Inf,Inf,k,J,as,bs,cs,skew,kurt,mt,st)$value^2


  #################################################
  #CTT Variance of Expected total scores
  #################################################
  Ex.theta=function(x,k,J,as,bs,cs){
    sapply(x,function(x) {
      val=numeric(length(x))
      for(j in 1:J){	for(m1 in 1:k){
        val=m1*Pij(x=x,m=m1,k,a=as[j],b=bs[j],c0=cs[m1],c1=cs[m1+1])+val
      }	}
      val
    })
  }

  Ex.theta.INT=function(x,k,J,as,bs,cs,skew,kurt,mt,st){
    Ex.theta(x,k,J,as,bs,cs)*dfleish(x,skew,kurt,mt,st)
  }
  Ex=integrate(Ex.theta.INT,-Inf,Inf,k,J,as,bs,cs,skew,kurt,mt,st)$value

  Var.x.INT=function(x,k,J,as,bs,cs,skew,kurt,mt,st,Ex){
    (Ex.theta(x,k,J,as,bs,cs)-Ex)^2*dfleish(x,skew,kurt,mt,st)
  }
  Vx=integrate(Var.x.INT,-Inf,Inf,k,J,as,bs,cs,skew,kurt,mt,st,Ex)$value

  rxx=Vx/(Vx+Vx.t)
  #################################################
  #Reliability
  #################################################
  c(rxx, thetaxx )
}
