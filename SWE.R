rm(list=ls(all=TRUE))

#-------------------------------------Parameters------------------------
nx<-50                  # grid size of x
ny<-31                  #grid size of y
n<-nx*ny
g<-9.8                  # gravitational constant
f<-1e-4
ka1<-1e+4
ka2<-5e+4
dt<-30                # hardwired timestep
L<-5e+5
D<-3e+5
dx<-L/nx
dy<-D/(ny-1)
nt<-8700 
nt_asm<-5820          # ͬ??48Сʱ??????ǰԤ??24h
dt_asm<-360           #ÿ3hͬ??һ??
nobs<-(nt_asm-60)%/%dt_asm+1    #??Ԥ??1h???ٿ?ʼͬ??        
ndim<-n*3
xobs<-ny*10*3                          # ֻ??h?Ĺ۲?
m<-100                                   # ??????????      
sigo<-c(sqrt(0.5),sqrt(0.5),1)                              # ?۲?????
sigi<-1                               # ??ʼ????
id<-c(5,9,16,21,25,30,31,33,44,49)       

#---------------------------------------------variables---------------------------------------
vf<-1
s<-matrix(,nrow=ndim,ncol=nt)        # ͬ????״̬??��????ֵ
s1<-matrix(,nrow=ndim,ncol=nt)       # ??ͬ??ʱ??Ԥ?????ϵľ?ֵ??ͬ??ʱ?̷????????ϵľ?ֵ
s2<-matrix(,nrow=ndim,ncol=nt)        
s3<-matrix(,nrow=ndim,ncol=nt)        
lamo<-array(NA,nobs)
y_adj<-array(NA,nobs)                    # ?Ż???????Ȼ???? 
y_0<-array(NA,nobs)                      # ?Ż?ǰ????Ȼ????
tobs<-rep(0,times=nt)
t_start<-60                             # ?й۲?????ʼʱ??
for(i in seq(t_start,nt_asm,by=dt_asm)) {tobs[i]<-1}     # tobs[i]Ϊ1ʱ??iʱ???й۲?
obs<-matrix(,nrow=xobs,ncol=nobs)        # ??ֵ????1???Ŷ????Ĺ۲?ֵ
obs_en<-array(NA,c(xobs,m,nobs))         # ?۲?ֵ????m???Ŷ????Ĺ۲?ֵensemble
do<-array(NA,c(xobs,m,nobs))
di1<-array(NA,c(xobs,m))                 # ?в?y-HX???Ż?ǰ??
di2<-array(NA,c(xobs,m))                 # ?в?y-HX???Ż?????       
s_en1<-array(NA,c(ndim,m))               # ?????Ŷ?????״̬??��ensemble???Ż?ǰ??
s_en2<-array(NA,c(ndim,m))               # ?????Ŷ?????״̬??��ensemble???Ż?????
s_a1<-array(NA,c(ndim,m))                # ??????ensemble???Ż?ǰ??
s_a2<-array(NA,c(ndim,m))                # ??????ensemble???Ż?????
s_m1<-matrix(,nrow=ndim,ncol=nobs)      # ͬ??ʱ?̱????????ϵľ?ֵ 
s_m2<-matrix(,nrow=ndim,ncol=nobs) 
a_m1<-matrix(,nrow=ndim,ncol=nobs)       # ͬ??ʱ?̷????????ϵľ?ֵ
a_m2<-matrix(,nrow=ndim,ncol=nobs)                       
z1<-matrix(,nrow=ndim,ncol=m) 
z2<-matrix(,nrow=ndim,ncol=m)      
#---------------------------------------------?۲?????(v??hȫ?۲?)------------------------------------------------
H<-matrix(0,nrow=xobs,ncol=ndim)         
for(i in 1:10)
  {jy<-(id[i]-1)*ny+1:ny
   jx<-((i-1)*ny+1):(i*ny)   
   for(ii in 1:ny) 
     {H[jx[ii],jy[ii]]<-1
      H[jx[ii]+ny*10,jy[ii]+n]<-1
      H[jx[ii]+ny*10*2,jy[ii]+n*2]<-1
     }
  } 
#-------------------------------------------------------END------------------------------------------------------
#---------------------------------------------------?۲?????????R------------------------------------------------
corr<-0.5
aa1<-matrix(,nrow=2,ncol=ny*10)         # ?۲????Ŀռ???ά????
for(i in 1:10)
  {iy<-((i-1)*ny+1):(i*ny)   
   for(ii in 1:ny) 
     {aa1[1,iy[ii]]<-ii
      aa1[2,iy[ii]]<-id[i]
     }
  } 
aa2<-diag(1,ny*10)
for(i in 1:(ny*10-1))
  {for(j in (i+1):(ny*10)) 
    {dd<-sqrt((aa1[1,i]-aa1[1,j])^2+(aa1[2,i]-aa1[2,j])^2)
     aa2[i,j]<-corr^dd
     aa2[j,i]<-aa2[i,j]
    }
  }
R<-matrix(0,nrow=xobs,ncol=xobs)
R[1:(ny*10),1:(ny*10)]<-aa2*vf*sigo[1]^2
R[(ny*10+1):(ny*20),(ny*10+1):(ny*20)]<-aa2*vf*sigo[2]^2
R[(ny*20+1):(ny*30),(ny*20+1):(ny*30)]<-aa2*vf*sigo[3]^2

bb1<-svd(R[1:(ny*10),1:(ny*10)]/vf)                             
sqrR_true1<-bb1$u%*%diag(sqrt(bb1$d),ny*10)            # ???ڲ????۲⣨??ʵ??R??
bb2<-svd(R[(ny*10+1):(ny*20),(ny*10+1):(ny*20)]/vf)                             
sqrR_true2<-bb2$u%*%diag(sqrt(bb2$d),ny*10)            # ???ڲ????۲⣨??ʵ??R??
bb3<-svd(R[(ny*20+1):(ny*30),(ny*20+1):(ny*30)]/vf)                             
sqrR_true3<-bb3$u%*%diag(sqrt(bb3$d),ny*10)            # ???ڲ????۲⣨??ʵ??R??
bb4<-svd(R)                               
sqrR<-solve(bb4$u%*%diag(sqrt(bb4$d),xobs))           # ???ڼ???????ʽ
inverR<-solve(R)
logR<-sum(log(bb4$d))
#------------------------------------------------------END-------------------------------------------------------

#-------------------------------------Generate true state-------------------------------------------------------
#--------------------Initialization----------------------------------
x<-array(rep(0:(nx-1),each=ny),dim<-c(ny,nx))/nx*L
y<-array(rep(0:(ny-1),times=nx),dim<-c(ny,nx))/(ny-1)*D

H0<-matrix(1,nrow=ny+2,ncol=nx+2)
U0<-matrix(0,nrow=ny+2,ncol=nx+2)
V0<-matrix(0,nrow=ny+2,ncol=nx+2)
Hx<-matrix(0,nrow=ny+2,ncol=nx+1)
Ux<-matrix(0,nrow=ny+2,ncol=nx+1)
Vx<-matrix(0,nrow=ny+2,ncol=nx+1)
Hy<-matrix(0,nrow=ny+1,ncol=nx+2)
Uy<-matrix(0,nrow=ny+1,ncol=nx+2)
Vy<-matrix(0,nrow=ny+1,ncol=nx+2)

h0<-50
h1<-5.5
h2<-3.325
H0[2:(ny+1),2:(nx+1)]<-h0+h1*tanh(9*(D/2-y)/(2*D))+h2/cosh(9*(D/2-y)/D)^2*sin(2*pi/L*x)
U0[2:(ny+1),2:(nx+1)]<--g/f*(-h1*(1-tanh(9*(D/2-y)/(2*D))^2)*9/(2*D)+2*h2/cosh(9*(D/2-y)/D)^2*tanh(9*(D/2-y)/D)*sin(2*pi/L*x)*9/D)
V0[2:(ny+1),2:(nx+1)]<-g/f*h2/cosh(9*(D/2-y)/D)^2*cos(2*pi/L*x)*2*pi/L

s[,1]<-c(as.vector(U0[2:(ny+1),2:(nx+1)]),as.vector(V0[2:(ny+1),2:(nx+1)]),as.vector(H0[2:(ny+1),2:(nx+1)]))

swe<-function(s_l,g_l,ka)
   {  # Reflective boundary conditions
	   H0[2:(ny+1),2:(nx+1)]<-array(s_l[(2*n+1):(3*n)],dim=c(ny,nx))
	   U0[2:(ny+1),2:(nx+1)]<-array(s_l[1:n],dim=c(ny,nx))*H0[2:(ny+1),2:(nx+1)]
	   V0[2:(ny+1),2:(nx+1)]<-array(s_l[(n+1):(2*n)],dim=c(ny,nx))*H0[2:(ny+1),2:(nx+1)]
	   
       H0[,1]<-H0[,nx+1]
       U0[,1]<-U0[,nx+1]
       V0[,1]<-V0[,nx+1]
       H0[,nx+2]<-H0[,2]
       U0[,nx+2]<-U0[,2]
       V0[,nx+2]<-V0[,2]
       H0[1,]<-H0[2,]
       U0[1,]<-U0[2,]
       V0[1,]<--V0[2,]
       H0[ny+2,]<-H0[ny+1,]
       U0[ny+2,]<-U0[ny+1,]
       V0[ny+2,]<--V0[ny+1,]

       # First half step
       # x direction
       i<-1:(ny+2)
       j<-1:(nx+1)
       # height
       Hx[i,j]<-(H0[i,j]+H0[i,j+1])/2-dt/(2*dx)*(U0[i,j+1]-U0[i,j])
       # x momentum
       Ux[i,j]<-(U0[i,j]+U0[i,j+1])/2-dt/(2*dx)*((U0[i,j+1]^2/H0[i,j+1]+g_l/2*H0[i,j+1]^2)-(U0[i,j]^2/H0[i,j]+g_l/2*H0[i,j]^2))
       # y momentum
       Vx[i,j]<-(V0[i,j]+V0[i,j+1])/2-dt/(2*dx)*((U0[i,j+1]*V0[i,j+1]/H0[i,j+1])-(U0[i,j]*V0[i,j]/H0[i,j]))
       
       # y direction
       i<-1:(ny+1)
       j<-1:(nx+2)
       # height
       Hy[i,j]<-(H0[i,j]+H0[i+1,j])/2-dt/(2*dy)*(V0[i+1,j]-V0[i,j])
       # x momentum
       Uy[i,j]<-(U0[i,j]+U0[i+1,j])/2-dt/(2*dy)*((V0[i+1,j]*U0[i+1,j]/H0[i+1,j])-(V0[i,j]*U0[i,j]/H0[i,j]))
       # y momentum
       Vy[i,j]<-(V0[i,j]+V0[i+1,j])/2-dt/(2*dy)*((V0[i+1,j]^2/H0[i+1,j]+g_l/2*H0[i+1,j]^2)-(V0[i,j]^2/H0[i,j]+g_l/2*H0[i,j]^2))
       
       # Second half step
       i<-2:(ny+1)
       j<-2:(nx+1)
       # height
       H0[i,j]<-H0[i,j]-(dt/dx)*(Ux[i,j]-Ux[i,j-1])-(dt/dy)*(Vy[i,j]-Vy[i-1,j])
       # x momentum
       U0[i,j]<-U0[i,j]+dt*f*V0[i,j]-(dt/dx)*((Ux[i,j]^2/Hx[i,j]+g_l/2*Hx[i,j]^2)-(Ux[i,j-1]^2/Hx[i,j-1]+g_l/2*Hx[i,j-1]^2))-(dt/dy)*((Vy[i,j]*Uy[i,j]/Hy[i,j])-(Vy[i-1,j]*Uy[i-1,j]/Hy[i-1,j]))
       # y momentum
       V0[i,j]<-V0[i,j]-dt*f*U0[i,j]-(dt/dx)*((Ux[i,j]*Vx[i,j]/Hx[i,j])-(Ux[i,j-1]*Vx[i,j-1]/Hx[i,j-1]))-(dt/dy)*((Vy[i,j]^2/Hy[i,j] + g_l/2*Hy[i,j]^2)-(Vy[i-1,j]^2/Hy[i-1,j]+g_l/2*Hy[i-1,j]^2))

       H0[,1]<-H0[,nx+1]
       U0[,1]<-U0[,nx+1]
       V0[,1]<-V0[,nx+1]
       H0[,nx+2]<-H0[,2]
       U0[,nx+2]<-U0[,2]
       V0[,nx+2]<-V0[,2]
       H0[1,]<-H0[2,]
       U0[1,]<-U0[2,]
       V0[1,]<--V0[2,]
       H0[ny+2,]<-H0[ny+1,]
       U0[ny+2,]<-U0[ny+1,]
       V0[ny+2,]<--V0[ny+1,]
	   
	   U0<-U0/H0
	   V0<-V0/H0
	   i<-2:(ny+1)
	   j<-2:(nx+1)
	   U0[i,j]<-U0[i,j]+dt/dx^2*ka*(U0[i,j+1]+U0[i,j-1]-2*U0[i,j])+dt/dy^2*ka*(U0[i+1,j]+U0[i-1,j]-2*U0[i,j])
	   V0[i,j]<-V0[i,j]+dt/dx^2*ka*(V0[i,j+1]+V0[i,j-1]-2*V0[i,j])+dt/dy^2*ka*(V0[i+1,j]+V0[i-1,j]-2*V0[i,j])
	   H0[i,j]<-H0[i,j]+dt/dx^2*ka*(H0[i,j+1]+H0[i,j-1]-2*H0[i,j])+dt/dy^2*ka*(H0[i+1,j]+H0[i-1,j]-2*H0[i,j])
      
       return(c(as.vector(U0[2:(ny+1),2:(nx+1)]),as.vector(V0[2:(ny+1),2:(nx+1)]),as.vector(H0[2:(ny+1),2:(nx+1)])))
    }
	
	for(t in 2:nt) {s[,t]<-swe(s[,t-1],g,ka1)}
#-----------------------------------------------------------end of generating true state------------------------------------------------

#---------------------------------------------------EnKF assimilation-----------------------------------------------------------
inverse<-function(zt)
  {
  inverR-inverR%*%H%*%zt%*%solve(diag(1,m)+t(H%*%zt)%*%inverR%*%(H%*%zt))%*%t(H%*%zt)%*%inverR
  }    

so<-H%*%s[,tobs==1]
ww<-matrix(,nrow=ny*10,ncol=nobs)
for(i in 1:nobs) {ww[,i]<-rnorm(ny*10,0,1)}
obs_err1<-sqrR_true1%*%ww
obs_err2<-sqrR_true2%*%ww
obs_err3<-sqrR_true3%*%ww
ww_en<-array(rnorm(ny*10*m*nobs,0,1),c(ny*10,m,nobs))*sqrt(vf) 
for(i in 1:nobs)
  {do[1:(ny*10),,i]<-sqrR_true1%*%ww_en[,,i]
   do[(ny*10+1):(ny*20),,i]<-sqrR_true2%*%ww_en[,,i]
   do[(ny*20+1):(ny*30),,i]<-sqrR_true3%*%ww_en[,,i]
  }
obs[1:(ny*10),]<-so[1:(ny*10),]+obs_err1
obs[(ny*10+1):(ny*20),]<-so[(ny*10+1):(ny*20),]+obs_err2
obs[(ny*20+1):(ny*30),]<-so[(ny*20+1):(ny*30),]+obs_err3                 
for(i in 1:m) {obs_en[,i,]<-obs+do[,i,]}                                          # ?۲?????m???Ŷ????γɹ۲?ensemble
ran<-matrix(rnorm(n*m,0,sigi),nrow=n)
for(j in 1:m) 
  {s_en1[(n*2+1):ndim,j]<-s[(n*2+1):ndim,1]+ran[,j]       # ???Ƴ?ֵ????m???Ŷ????γɳ?ֵensemble
	   H0[2:(ny+1),2:(nx+1)]<-array(s_en1[(2*n+1):(3*n),j],dim=c(ny,nx))
	 H0[,1]<-H0[,nx+1]
	 H0[,nx+2]<-H0[,2]
	 H0[1,]<-H0[2,]
	 H0[ny+2,]<-H0[ny+1,]
	 k<-2:(ny+1)
	 l<-2:(nx+1)
     U0[k,l]<--g/f*(H0[k+1,l]-H0[k-1,l])/(2*dy)
     V0[k,l]<-g/f*(H0[k,l+1]-H0[k,l-1])/(2*dx)
     s_en1[1:(2*n),j]<-c(as.vector(U0[2:(ny+1),2:(nx+1)]),as.vector(V0[2:(ny+1),2:(nx+1)]))
  }
s_en2<-s_en1                                                           
s1[,1]<-apply(s_en1,1,mean)     
s2[,1]<-apply(s_en2,1,mean)
s3[,1]<-s[,1]  

for(t in 1:(nt-1)) {s3[,t+1]<-swe(s3[,t],g,ka2)}


num<-0                                                           
for(t in 1:(nt-1)) {
  print(t)
 if(tobs[t]!=1)                                         # ??????ʱ???޹۲⣬ֻԤ????ͬ??
  {for(j in 1:m)
    {s_en1[,j]<-swe(s_en1[,j],g,ka2)
     s_en2[,j]<-swe(s_en2[,j],g,ka2)      
    }
   s1[,t+1]<-apply(s_en1,1,mean) 
   s2[,t+1]<-apply(s_en2,1,mean)     
  } 
  else                                                  # ??????ʱ???й۲⣬ͬ??
     {num<-num+1          
      di1<-obs_en[,,num]-H%*%s_en1       
      di2<-obs_en[,,num]-H%*%s_en2
      s_m1[,num]<-apply(s_en1,1,mean)
      s_m2[,num]<-apply(s_en2,1,mean)
      di_m1<-obs[,num]-H%*%s_m1[,num]     
      di_m2<-obs[,num]-H%*%s_m2[,num]      
      for(j in 1:m) {z1[,j]<-(s_en1[,j]-s_m1[,num])*sqrt(1/(m-1))}
      for(j in 1:m) {z2[,j]<-(s_en2[,j]-s_m2[,num])*sqrt(1/(m-1))}
      P1<-z1%*%t(z1)
      P2<-z2%*%t(z2)
#--------------------------------------------------------------------------δ?Ż?ǰ?ĵ?ͬ??????---------------------------------------------------
      KM1<-P1%*%t(H)%*%inverse(z1)                      # ????????K
      s_a1<-s_en1+KM1%*%di1                                                   # ??????
      for(j in 1:m) {s_en1[,j]<-swe(s_a1[,j],g,ka2)}
      mul<-1
      DHZ<-array(0,xobs)
      DHZ[1:m]<-svd(sqrR%*%H%*%z1)$d
      for(i in 1:xobs) {mul<-mul*(1+DHZ[i]^2)}      
      y_0[num]<-logR+log(mul)+t(di_m1)%*%inverse(z1)%*%di_m1
      a_m1[,num]<-apply(s_a1,1,mean)              # ????????ֵ                               
      s1[,t]<-apply(s_a1,1,mean)
      s1[,t+1]<-apply(s_en1,1,mean)
#-----------------------------------------------------------------------------end----------------------------------------------------------------
#--------------------------------------------------------------------------?Ż?????ͬ??????-------------------------------------------------------         
      Po<-P2
      di_mo<-di_m2 
      DHZo<-array(0,xobs)
      DHZo[1:m]<-svd(sqrR%*%H%*%z2)$d   
      likelio<-function(xmin)                                        # ??Ȼ????
         {y1<-xmin
          mulo<-1
          for(i in 1:xobs) {mulo<-mulo*(1+(y1*DHZo[i])^2)}      
          logR+log(mulo)+t(di_mo)%*%inverse(z2*y1)%*%di_mo
         }      
      lamo[num]<-optimize(likelio,c(0.1,20),tol=0.01)$minimum
      y_adj[num]<-optimize(likelio,c(0.1,20),tol=0.01)$objective
      P_f<-lamo[num]^2*P2                                # ?Ż?????Ԥ??????Э????????     
      s_p<-matrix(,nrow=ndim,ncol=m)
      for(j in 1:m) {s_p[,j]<-s_m2[,num]+lamo[num]*(s_en2[,j]-s_m2[,num])}
      di2<-obs_en[,,num]-H%*%s_p
      KM2<-P_f%*%t(H)%*%inverse(z2*lamo[num])     # ????????K
      s_a2<-s_p+KM2%*%di2                                # ??????
      s_en2<-s_p 
      for(j in 1:m) {s_en2[,j]<-swe(s_a2[,j],g,ka2)}
      a_m2[,num]<-apply(s_a2,1,mean)                           # ????????ֵ 
      s2[,t]<-apply(s_a2,1,mean)
      s2[,t+1]<-apply(s_en2,1,mean)                           
#-----------------------------------------------------------------------------end----------------------------------------------------------------
     write(lamo[num],file="Output/lamo.txt",append=TRUE,sep="\t") 
   } 
 }
#--------------------------------------------------------------end of EnKF assimilation-----------------------------------------------------------
save.image("Output/930obs.RData")

plot(apply(s - s1, 2, function(x){sqrt(mean(x ^ 2))}), type = 'l')
lines(apply(s - s2, 2, function(x){sqrt(mean(x ^ 2))}), col = 'red')
