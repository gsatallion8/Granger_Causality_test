library(complexplus)

Vhat <- function(Data,Order)
  {
  # Function to obtain V_ij matrix (P*P) from the estimate of covariance matrix (M^2P*M^2P) of
  # VAR model coefficients
  # Reference: Michael Eichler (2005) and Lutkepohl H (1993)
  # Usage: [Vmat] = Vhat(Data, Order)
  # Sudhakar Kathari, Srikanth Sarma, November 11, 2017.

  # Read input arguments
  Z = Data
  P = Order
  
  # Size of the data
  Size = dim(Z)
  N = Size[1]
  M = Size[2]
  
  # The regressor matrix
  
   BigZ = {}
   BigZ1 = {}
   
   for (k in 1:P){
     BigZ = cbind(BigZ, Z[k:(N-P+k-1),])
   }
  
   # Re-arrange the regressor matrix
   for (j in 1:M){
     BigZ1 = cbind(BigZ1,BigZ[,(seq(j,M*P,M))])
   }
   # The Y matrix   
   Y = Z[(P+1):N,]
   
   # Regressors covariance matrix and its inverse
   SigmaR = (t(BigZ1)%*%BigZ1)/N
   invSigmaR = qr.solve(SigmaR)
   
   # Estimate VAR model coefficient matrices
   Ahat = (qr.solve(t(BigZ1)%*%BigZ1)%*%t(BigZ1))%*%Y
   
   # Innovations covariance matrix
   SigmaE = (1/(N-P))*t((Y-BigZ1%*%Ahat))%*%(Y-BigZ1%*%Ahat)
   
   # Covariance matrix of VAR model coefficients
   SigmaA = kronecker(invSigmaR,SigmaE)
   
   
   # Estimation of V matrix
  Vmatrix = array(0,dim = c(P,P,M,M))
   for (i in seq(1,M)){
     for(j in seq(1,M))
     Vmatrix[,,i,j] = SigmaA[(((i-1)*M+j-1)*P+1):(((i-1)*M+j)*P),(((i-1)*M+j-1)*P+1):(((i-1)*M+j)*P)]
   }
  A<-array(0,dim=c(M,M,P))
  for(i in seq(1,M))
  {
    for(j in seq(1,M))
    {
      for(k in seq(1,P))
      {
        A[i,j,k]<-Ahat[P*(j-1)+k,i]
      }
    }
  }
  c<-list()
  c[[1]]<-A
  c[[2]]<-Vmatrix
   return(c)
}


GrangerTest_timeDomain<-function(Zk,P)
{
  c<-Vhat(Zk,P)
  A<-c[[1]]
  V<-c[[2]]
  size=dim(Zk)
  N=size[1]
  M=size[2]
  
  S<-array(0,dim=c(M,M))
  for(i in seq(1,M))
  {
    for(j in seq(1,M))
    {
      S[i,j]=N*t(A[i,j,])%*%qr.solve(V[,,i,j])%*%A[i,j,]
    }
  }
  return(S)
}


FRF_Matrix<-function(A,w)
{
  P=dim(A)[3]
  M=dim(A)[1]
  Abar<-diag(M)
  for(k in seq(1,P))
  {
    Abar<-Abar-A[,,k]*exp(-1i*k*w)
  }
  
  return(Abar)
}


coherency_matrix<-function(Abar)
{
  size=dim(Abar)
  M=size[1]
  
  Psi<-array(0,dim=c(M,M))
  h<-array(0,dim=c(M,M))
  sum<-array(0,dim=M)
  detA<-Det(Abar)
  
  if(M>3)
  {
    for(i in seq(1,M))
    {
      for(j in seq(1,M))
      {
        if(i==j)
        {
          h[i,j]=Det(Abar[-i,-j])/detA
        }
        else
        {
          h[i,j]=-Abar[i,j]*Det(Abar[-c(i,j),-c(i,j)])/detA
        }
      }
    }
    for(j in seq(1,M))
    {
      for(i in seq(1,M))
      {
        sum[j]=sum[j]+Mod(h[i,j])*Mod(h[i,j])
      }
      sum[j]=sqrt(sum[j])
    }
    for(i in seq(1,M))
    {
      for(j in seq(1,M))
      {
        Psi[i,j]=h[i,j]/sum[i]
      }
    }
  }
  else
  {
    for(i in seq(1,M))
    {
      for(j in seq(1,M))
      {
        if(i==j)
        {
          h[i,j]=Det(Abar[-i,-j])/detA
        }
        else
        {
          h[i,j]=-Abar[i,j]*Abar[-c(i,j),-c(i,j)]/detA
        }
      }
    }
    for(j in seq(1,M))
    {
      for(i in seq(1,M))
      {
        sum[j]=sum[j]+Mod(h[i,j])*Mod(h[i,j])
      }
      sum[j]=sqrt(sum[j])
    }
    
    for(i in seq(1,M))
    {
      for(j in seq(1,M))
      {
        Psi[i,j]=h[i,j]/sum[j]
      }
    }
  }

  
  return(Psi)
}


DPF<-function(Zk,P,Nf,plot)
{
  #Nf=N for general fft
  
  size=dim(Zk)
  N=size[1]
  M=size[2]
  
  Psi<-array(0,dim=c(M,M,Nf))
  wseq<-seq(0,2*pi,2*pi/(Nf-1))
  
  c<-Vhat(Zk,P)
  A<-c[[1]]
  V<-c[[2]]
    
  k=1
  for(w in wseq)
  {
    Abar<-FRF_Matrix(A,w)
    Psi[,,k]<-coherency_matrix(Abar)
    k=k+1
  }
  
  #Plotting
  if(plot==TRUE)
  {
    par(mfrow=c(M,M))
    for(i in seq(1,M))
    {
      for(j in seq(1,M))
      {
        plot(wseq,Mod(Psi[i,j,])^2,type='l',ylab = expression(Psi^2),xlab = expression(omega))
      }
    }
  }
  
  return(Psi)
}


Surrogate_Generation<-function(Zk,R,P,Nf)
{
  Zw<-Mod(mvfft(Zk))
  size=dim(Zk)
  N=size[1]
  M=size[2]
  Zk_surr<-array(0,dim=c(N,M,R))
  Zw_surr<-array(0,dim=dim(Zw))
  Psi_surr<-array(0,dim=c(M,M,Nf,R))
  for(i in seq(1,R))
  {
    for(y in seq(1,M))
    {
      for(x in seq(1,N))
      {
        Zw_surr[x,y]=Zw[x,y]*exp(1i*runif(1,min=0,max=2*pi))
      }
      Zk_surr[,y,i]<-Re(fft(Zw_surr[,y],inverse=TRUE))/length(Zw_surr[,y])
    }
    Psi_surr[,,,i]<-DPF(Zk_surr[,,i],P,Nf,plot = FALSE)
  }
  
  return(Psi_surr)
}


A1<-array(0,dim=c(3,3))
A2<-array(0,dim=c(3,3))
A1[1,1]=0.3
A1[2,1]=0.6
A1[2,2]=0.4
A1[3,2]=0.4
A1[3,3]=0.5
A2[1,1]=0.2
A2[2,2]=0.5
A2[3,2]=0.3
A2[3,3]=0.4
sigma<-diag(3)

z<-VARMAsim(2000,arlags = 2,phi=cbind(A1,A2),sigma = sigma)

z<-z$series

model<-VAR(z,lag.max = 10,ic='SC')
P<-model$p


GrangerTest_timeDomain(z,P)<qchisq(0.95,df=P)

Psi<-DPF(z,3,1000,plot = TRUE)

R=100
Nf=100
Psi_surr<-Surrogate_Generation(z,R,P,Nf)

size=dim(Psi_surr)
M=size[1]

Psi<-DPF(z,3,1000,plot = FALSE)

Psi_95<-array(0,dim=c(M,M,Nf))
non_Cause<-array(TRUE,dim=c(M,M))

for(i in seq(1,M))
{
  for(j in seq(1,M))
  {
    for(k in seq(1,Nf))
    {
      Psi_95[i,j,k]<-quantile(Mod(Psi_surr[i,j,k,])^2,probs=0.95)
      if(Mod(Psi[i,j,k])^2>Psi_95[i,j,k])
      {
        non_Cause[i,j]=FALSE
        break
      }
    }
  }
}

non_Cause

model<-VAR(Zk,lag.max = 10,ic='SC')
P<-model$p
GrangerTest_timeDomain(Zk,P+1)<qchisq(0.95,df=P+1)

Psi<-DPF(Zk,P+1,1000,plot = TRUE)

R=100
Nf=100
Psi_surr<-Surrogate_Generation(Zk,R,P+1,Nf)

size=dim(Psi_surr)
M=size[1]

Psi_95<-array(0,dim=c(M,M,Nf))
non_Cause<-array(TRUE,dim=c(M,M))

for(i in seq(1,M))
{
  for(j in seq(1,M))
  {
    for(k in seq(1,Nf))
    {
      Psi_95[i,j,k]<-quantile(Mod(Psi_surr[i,j,k,])^2,probs = 0.95)
      if(Mod(Psi[i,j,k])^2>Psi_95[i,j,k])
      {
        non_Cause[i,j]=FALSE
        break
      }
    }
  }
}

non_Cause