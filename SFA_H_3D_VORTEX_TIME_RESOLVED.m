%SFA_H_3D_VORTEX_TIME_RESOLVED.m: SFA, real hydrogen atom
%for Calvin He on Aug.11 2017
clear,clc,tic
I1=sqrt(-1);

%Two oppositely circularly polarized time-delayed as pulses: PRL 115 113004 2015
nt=1000;
W=2/27.2;    
T=2*pi/W;  
tau=1*T;
E0=sqrt(3.5e14/3.5e16); 
Ip=13.6/27.2;  
ksi=1;
dt=8*tau/(nt-1);   
t(nt)=NaN;
td(nt)=NaN;
Ex(nt)=NaN;
Ey(nt)=NaN;
F(nt)=NaN;
Fd(nt)=NaN;
for i=1:nt
t(i) =(i-3.5*nt/10)*dt;
td(i)=(i-6.5*nt/10)*dt;
F(i)=exp(-t(i)^2/tau^2);
Fd(i)=exp(-td(i)^2/tau^2);
Ex(i)=   F(i)*E0    /(1+ksi^2)^0.5*cos(1.0*W*t(i)+0) + 1*Fd(i)*E0    /(1+ksi^2)^0.5*cos(1.0*W*td(i)+0*pi/2);  
Ey(i)=   F(i)*E0*ksi/(1+ksi^2)^0.5*sin(1.0*W*t(i)+0) - 1*Fd(i)*E0*ksi/(1+ksi^2)^0.5*sin(1.0*W*td(i)+0*pi/2);  
end

%laser-pulse vector potential
Ax(nt)=NaN;
Ay(nt)=NaN;
for i=1:nt
jkx=0;
jky=0;
for j=1:i
jkx=jkx-Ex(j)*dt;
jky=jky-Ey(j)*dt;
end
Ax(i)=jkx;
Ay(i)=jky;
end

%Momentum ranges
nx=100;              
Px(nx)=NaN;
dpx=3/(nx-1); 
for i=1:nx
Px(i)=(i-nx/2-1/2)*dpx;    
end
ny=100;
Py(ny)=NaN;
dpy=3/(ny-1);
for i=1:ny
Py(i)=(i-ny/2-1/2)*dpy;    
end

%time-resolved momentum spectra
Mx(nx,ny)=NaN;
My(nx,ny)=NaN;
P(nt,nx,ny)=NaN; %This function is the time-resolved momenta, #frames=nt/50+1. 

tic
for nt0=1:50:nt
for ix=1:nx   
for iy=1:ny    
    
    Mx(ix,iy)=(sqrt(-1)*2^(3.5)*(2*Ip)^(5/4)/pi)*Px(ix)/(Px(ix)^2+Py(iy)^2+2*Ip)^3;
    My(ix,iy)=(sqrt(-1)*2^(3.5)*(2*Ip)^(5/4)/pi)*Py(iy)/(Px(ix)^2+Py(iy)^2+2*Ip)^3; 
    jkM=0;    
    for it=1:nt0        
        phase=0;
        for it2=it:nt0
        phase=phase+dt*( (Px(ix)+Ax(it2))^2/2+(Py(iy)+Ay(it2))^2/2+Ip );
        end    
        jkM=jkM+sqrt(-1)*dt*Ex(it)*Mx(ix,iy)*exp(-sqrt(-1)*phase);
        jkM=jkM+sqrt(-1)*dt*Ey(it)*My(ix,iy)*exp(-sqrt(-1)*phase);
    end
    P(nt0,ix,iy)=jkM;
end %iy
end %ix
end %nt0
toc/60
%end--------------------------
