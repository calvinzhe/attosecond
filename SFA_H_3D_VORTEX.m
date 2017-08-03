%SFA_H_3D_VORTEX.m: SFA, real hydrogen atom, Feb20 2017
%works
clear,clc,tic
close all
I1=sqrt(-1);
tic
%Two oppositely circularly polarized time-delayed as pulses: PRL 115 113004 2015
nt=100;
W=3/27.2;    %1.323
T=12*pi/W;  %
%T=2000; %pulse width in atomic units
tau=T;
E0=sqrt(3.5e14/3.5e16);
Ip=12.13/27.2;  %He=0.9037 for He
ksi=1;
dt=8*tau/(nt-1);   %8*
t(nt)=NaN;
td(nt)=NaN;
Ex(nt)=NaN;
Ey(nt)=NaN;
F(nt)=NaN;
Fd(nt)=NaN;
parfor i=1:nt
    t(i) =(i-3.0*nt/10)*dt;
    td(i)=(i-7.0*nt/10)*dt;
    F(i)=exp(-t(i)^2/T^2);
    Fd(i)=exp(-td(i)^2/T^2);
    Ex(i)=   F(i)*E0    /(1+ksi^2)^0.5*cos(-1.0*W*t(i)) + 1*Fd(i)*E0    /(1+ksi^2)^0.5*cos(1.0*W*td(i)+0*pi/2);
    Ey(i)=   F(i)*E0*ksi/(1+ksi^2)^0.5*sin(-1.0*W*t(i)) + 1*Fd(i)*E0*ksi/(1+ksi^2)^0.5*sin(1.0*W*td(i)+0*pi/2);
end

figure, plot(t,Ex);

%Vector potential
Ax(nt)=NaN;
Ay(nt)=NaN;
parfor i=1:nt
    jkx=0;
    jky=0;
    for j=1:i
        jkx=jkx-Ex(j)*dt;
        jky=jky-Ey(j)*dt;
    end
    Ax(i)=jkx;
    Ay(i)=jky;
end

%Momentum
nx=100;
Px(nx)=NaN;
dpx=3/(nx-1);
parfor i=1:nx
    Px(i)=(i-nx/2-1/2)*dpx;
end
ny=100;
Py(ny)=NaN;
dpy=3/(ny-1);
parfor i=1:ny
    Py(i)=(i-ny/2-1/2)*dpy;
end

%Momentum spectra
Mx(nx,ny)=NaN;
My(nx,ny)=NaN;
P(nx,ny)=NaN;
tic
for it=1:nt
    parfor ix=1:nx
        for iy=1:ny

            %K=sqrt(Px(ix)^2+Py(iy)^2);
            %Mx(ix,iy)=sqrt(128/3)*sqrt(K)*((K-I1)/(-K-I1))^(I1/k)*sqrt(1+coth(pi/K));
            %
            Mx(ix,iy)=(sqrt(-1)*2^(3.5)*(2*Ip)^(5/4)/pi)*Px(ix)/(Px(ix)^2+Py(iy)^2+2*Ip)^3;
            My(ix,iy)=(sqrt(-1)*2^(3.5)*(2*Ip)^(5/4)/pi)*Py(iy)/(Px(ix)^2+Py(iy)^2+2*Ip)^3;
            jkM=0;
            phase=0;
            for it2=it:nt
                %phase=phase+dt*( (Px(ix)-Ax(it2))^2/2+(Py(iy)-Ay(it2))^2/2+Ip +0.1*(F(it)+Fd(it)) );
                phase=phase+dt*( (Px(ix)+Ax(it2))^2/2+(Py(iy)+Ay(it2))^2/2+Ip );
            end
            jkM=jkM+sqrt(-1)*dt*Ex(it)*Mx(ix,iy)*exp(-sqrt(-1)*phase);
            jkM=jkM+sqrt(-1)*dt*Ey(it)*My(ix,iy)*exp(-sqrt(-1)*phase);
            P(ix,iy)=P(ix,iy)+jkM;
        end
    end
    it
    fig = figure;
    imagesc(Px,Py,abs(P').^2);
    set(gca,'YDir','normal');
    xlabel('Px (a.u.)');
    ylabel('Py (a.u.)');
    axis([-1 1 -1 1]);
    title('Photoelectron Momentum Distribution');
    saveas(fig,strcat('./CCP_Time_Animation/',num2str(it),'.png'));
end


toc
% figure
% pcolor(Px,Py,abs(P').^2);
% shading interp;
% xlabel('Px (a.u.)');
% ylabel('Py (a.u.)');

% %ADK probability amplitudes at field extrema
% omega=9.1e-31*(1.6e-19)^4/(4*pi*8.85e-12)^2/(1.055e-34)^3;
% Ea=(9.1e-31)^2*(1.6e-19)^5/(4*pi*8.85e-12)^3/(1.055e-34)^4;
% rH=13.6/13.6;
% alpha=4*omega*(rH)^2.5*Ea;
% beta=2*(rH)^1.5*Ea/3;
%
% PES(nz,nx)=NaN;
% for iz=1:nz
% for ix=1:nx
%
%     jpes=0;
%     for it=1:nt
%     Wt=alpha/(abs(Ea*EZ(it))+1000)*exp(-beta/(abs(Ea*EZ(it))+1000));
%     St=0;
%     JW=0;
%     for j=it+1:nt
%     St=St+dt*((PZ(iz)-AZ(it))^2/2+(PX(ix)-AX(it))^2/2+Ip);
%     JW=JW+exp(-abs(Wt)^2);
%     end
%     jpes=jpes+JW*Wt*exp(-sqrt(-1)*St);
%     end
%     PES(iz,ix)=jpes;
% end
% iz
% end
% toc/60
% PES=PES/max(max(abs(PES)));
% contour(PX,PZ,(1-abs(PES')).^2,100);
% contour(PX,PZ,abs(PES').^2,100);
% xlabel('Pz (a.u.)');
% ylabel('Px (a.u.)');
%end--------------------------
