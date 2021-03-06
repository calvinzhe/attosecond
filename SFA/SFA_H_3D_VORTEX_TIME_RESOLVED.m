%SFA_H_3D_VORTEX_TIME_RESOLVED.m: SFA, real hydrogen atom
%for Calvin He on Aug.11 2017
clear,clc,tic
close all
I1=sqrt(-1);

tic

%Two oppositely circularly polarized time-delayed as pulses: PRL 115 113004 2015
nt=500;
W=0.0577;
lambda = (3*10^8)/(W/(2*pi))*(2.419*10^-17)/10^-9;
nc=6;
T=nc*pi/W;  
tau=3*T;
phi1=0;
phi2=0;
E0=sqrt(3.5e15/3.5e16); 
Ip=13.6/27.2;  
ksi=1;
dt=T/nt;   
t(nt)=NaN;
td(nt)=NaN;
Ex(nt)=NaN;
Ey(nt)=NaN;
F(nt)=NaN;
Fd(nt)=NaN;
dir = './Lin_Ident_time_resolved/';

polarization = 'Linear XX'; %Polarization variable (see below)
p1x = 1; p2x = 1; p1y = 1; p2y = 1; p1w = 1; p2w = 1; %Circular RR Default
if strcmp(polarization, 'Circular LR')
    p1w = -1; p2w = 1;
elseif strcmp(polarization, 'Circular RL')
    p1w = 1; p2w = -1;
elseif strcmp(polarization, 'Linear XX')
    p1y = 0; p2y = 0;
elseif strcmp(polarization, 'Linear YY')
    p1x = 0; p2x = 0;
elseif strcmp(polarization, 'Linear XY')
    p1y = 0; p2x = 0;
elseif strcmp(polarization, 'Linear YX')
    p1x = 0; p2y = 0;
elseif strcmp(polarization, 'Circular LL')
    p1w = -1; p2w = -1;
elseif strcmp(polarization, 'Circular RR')
    p1w = 1; p2w = 1;
end

for i=1:nt
    %
    t(i) = 4*(2*i-nt-1)*dt;
    F(i)=exp(-(t(i)+tau/2)^2/T^2);
    Fd(i)=exp(-(t(i)-tau/2)^2/T^2);
    %}
    %{
    t(i) =2*(2*i-nt-1)*dt;
    if (t(i)+tau/2 > T/2) || (t(i)+tau/2 < -T/2)
        F(i) = 0;
    else
        F(i)= cos(pi*(t(i)+tau/2)/T)^2;
    end %Envelop function for 1st pulse
    if (t(i)-tau/2 > T/2) || (t(i)-tau/2 < -T/2)
        Fd(i) = 0;
    else
        Fd(i)= cos(pi*(t(i)-tau/2)/T)^2;
    end %Envelop fuction for the 2nd (delayed) pulse
    %}
    Ex(i)=   p1x*F(i)*E0    /(1+ksi^2)^0.5*cos(1.0*p1w*W*(t(i)+tau/2)+phi1) + p2x*Fd(i)*E0    /(1+ksi^2)^0.5*cos(1.0*p2w*W*(t(i)-tau/2)+phi2);
    Ey(i)=   p1y*F(i)*E0*ksi/(1+ksi^2)^0.5*sin(1.0*p1w*W*(t(i)+tau/2)+phi1) + p2y*Fd(i)*E0*ksi/(1+ksi^2)^0.5*sin(1.0*p2w*W*(t(i)-tau/2)+phi2);
end

%%

plotrange = linspace(min(t),max(t),nt);

figA = figure;
plot(plotrange,E0*F/sqrt(2),'b', ...
        plotrange,E0*Fd/sqrt(2),'r', ...
        plotrange,sqrt(Ex.^2+Ey.^2),'g');
legend('1st Pulse Envelope','2nd Pulse Envelope','E-field Magnitude');
xlabel('time');
axis([min(plotrange) max(plotrange) 0 1]);

figB = figure;
plot3(plotrange,Ex,Ey);
xlabel('time');
ylabel('Ex');
zlabel('Ey');
set(gca, 'Ydir','reverse');
axis([min(plotrange) max(plotrange) -1 1 -1 1]);

saveas(figA, strcat(dir,'Envelope_E-field_Magnitude.png'));
saveas(figB, strcat(dir,'E-field.png'));

%%
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
nx=200;              
Px(nx)=NaN;
dpx=4/(nx-1); 
for i=1:nx
    Px(i)=(i-nx/2-1/2)*dpx;    
end
ny=200;
Py(ny)=NaN;
dpy=4/(ny-1);
for i=1:ny
    Py(i)=(i-ny/2-1/2)*dpy;    
end

%time-resolved momentum spectra
Mx(nx,ny)=NaN;
My(nx,ny)=NaN;
P(nt,nx,ny)=NaN; %This function is the time-resolved momenta, #frames=nt/50+1. 

tic
for nt0=1:nt
    parfor ix=1:nx   
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
    nt0
    save(strcat(dir,'data.mat'),'P','nx','ny','Px','Py','Ex','Ey','plotrange','tau','polarization','nt','phi1','phi2','lambda','T','nc');
end %nt0

%%
parfor nt0=1:nt
    if nt0 < 10
        zeroes = '000';
    elseif nt0 < 100
        zeroes = '00';
    elseif nt0 < 1000
        zeroes = '0';
    else
        zeroes = '';
    end
    PAmp = reshape(P(nt0,:,:),[nx,ny]);
    PAmp = abs(PAmp').^2;
    
    fig = figure('visible','on');
    
    ax1 = axes('Position',[0.1 0.1 .85 .85]);
    imagesc(ax1,Px,Py,PAmp);
    axis([-2 2 -2 2]);
    set(gca, 'YDir','normal');
    colorbar;
    caxis([0 5]);
    xlabel('Px (a.u.)');
    ylabel('Py (a.u.)');
    title('Photoelectron Momentum Distribution');
    str1 = {strcat('$$\tau = ', num2str(round(tau,1)), '\ au$$'), ...
            strcat('$$T = ', num2str(round(T,1)), '\ au$$'), ...
            strcat('$$\lambda = ', num2str(round(lambda,1)),'\ nm$$')};
    text(-1.95,1.65,str1,'Interpreter','latex','BackgroundColor','yellow');
    str2 = {strcat('$$\phi_1 = ', num2str(round(phi1/pi,1)), '\ \pi$$'), ...
            strcat('$$\phi_2 = ', num2str(round(phi2/pi,1)), '\ \pi$$')};
    text(-1.95,-1.7,str2,'Interpreter','latex','BackgroundColor','yellow');
    str3 = {strcat('$$n_c = ', num2str(nc),'$$')};
    text(1.5,1.85,str3,'Interpreter','latex','BackgroundColor','yellow');
    str4 = {strcat('$$Polarization = ', polarization,'$$')};
    text(0.25,-1.85,str4,'Interpreter','latex','BackgroundColor','yellow');
    pbaspect([1 1 1]);
    
    ax2 = axes('Position',[0.625 0.175 0.15 0.15]);
    plot3(ax2,plotrange(nt0:nt),Ex(nt0:nt),Ey(nt0:nt));
    set(gca, 'Ydir','reverse','xtick',[],'ytick',[],'ztick',[]);
    axis([plotrange(nt0) plotrange(nt0)+max(plotrange)-min(plotrange) -.75 .75 -.75 .75]);
    
    saveas(fig,strcat(dir,'Momentum_Distribution_w_E-field_',zeroes,num2str(nt0),'.png'));
    nt0
end

toc
%end--------------------------

