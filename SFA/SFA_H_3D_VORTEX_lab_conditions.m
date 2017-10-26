%SFA_H_3D_VORTEX_cos2_envelop.m: SFA, real hydrogen atom, Feb20 2017
%works
clear,clc,tic
close all
I1=sqrt(-1);
opengl('save', 'software');

tic

%Two oppositely circularly polarized time-delayed as pulses: PRL 115 113004 2015
nt=500;
W= 0.0577; %3/27.2;    %0.1103 angular frequency => Period = 2*pi/W = 56.9675 => F = 0.0176
lambda = (3*10^8)/(W/(2*pi))*(2.419*10^-17)/10^-9;  %wavelength in nm
nc = 6;
T=nc*pi/W;  %Pulse width = nc*(optical period) = nc*(pi/W) for Gaussian
%T=2000; %pulse width in atomic units
tau=4*T;       %Pulse separation
phi1 = 0;
phi2 = 0;
E0=sqrt(3.5e14/3.5e16);
Ip=12.13/27.2;  %Ionization potential in au (1 Hartree = 27.2eV)  |  He=0.9037 for He
ksi=1;
dt=T/nt;   %Time increment for loop iteration
t(nt)=NaN;
td(nt)=NaN;
Ex(nt)=NaN;
Ey(nt)=NaN;
F(nt)=NaN;
Fd(nt)=NaN;
dir = './Lab_Conditions_Lin_XX_phi_sweep_2/';

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

n_phi = 150;
phi_list = zeros(n_phi,1);
for i_phi=100:100
    phi_list(i_phi) = (i_phi-1)*2*pi/n_phi-pi;
    phi2 = phi_list(i_phi);
    parfor i=1:nt
        t(i) = 4*(2*i-nt-1)*dt;
        F(i)=exp(-(t(i)+tau/2)^2/T^2);
        Fd(i)=exp(-(t(i)-tau/2)^2/T^2);
        Ex(i)=   p1x*F(i)*E0    /(1+ksi^2)^0.5*cos(1.0*p1w*W*(t(i)+tau/2)+phi1) + p2x*Fd(i)*E0    /(1+ksi^2)^0.5*cos(1.0*p2w*W*(t(i)-tau/2)+phi2);
        Ey(i)=   p1y*F(i)*E0*ksi/(1+ksi^2)^0.5*sin(1.0*p1w*W*(t(i)+tau/2)+phi1) + p2y*Fd(i)*E0*ksi/(1+ksi^2)^0.5*sin(1.0*p2w*W*(t(i)-tau/2)+phi2);
    end

    plotrange = linspace(min(t),max(t),nt);
    
    if i_phi < 10, zeroes = '00';
    elseif i_phi < 100, zeroes = '0';
    elseif i_phi < 1000, zeroes = '';
    end

    %{
    figA = figure('visible','off');
    plot(plotrange,E0*F/sqrt(2),'b', ...
            plotrange,E0*Fd/sqrt(2),'r', ...
            plotrange,sqrt(Ex.^2+Ey.^2),'g');
    legend('1st Pulse Envelope','2nd Pulse Envelope','E-field Magnitude');
    xlabel('time');
    axis([min(plotrange) max(plotrange) 0 1]);

    figB = figure('visible','off');
    plot3(plotrange,Ex,Ey);
    xlabel('time');
    ylabel('Ex');
    zlabel('Ey');
    axis([min(plotrange) max(plotrange) -1 1 -1 1]);
    
    saveas(figA, strcat(dir,'Envelope_E-field_Magnitude',zeroes,num2str(i_phi),'.png'));
    saveas(figB, strcat(dir,'E-field_',zeroes,num2str(i_phi),'.png'));
    %}
    

    %Vector potential
    % E(t) = -dA(t)/dt  =>  A(t) = -Integral (t0 to t) E(t')dt'
    Ax(nt)=NaN;
    Ay(nt)=NaN;
    parfor i=1:nt
        jkx=0;
        jky=0;
        for j=1:i
            jkx=jkx-Ex(j)*dt;   %Subtract Ex(t')dt' for t'=j*dt from jkx (adding one slice of the integral)
            jky=jky-Ey(j)*dt;   %Do the same for jky with Ey
        end
        Ax(i)=jkx;  %Set Ax(t) = jkx = - Integral (1 to t) Ex(t')dt'
        Ay(i)=jky;  %Do the same for jky with Ey
    end

    %Momentum
    nx=200;
    Px(nx)=NaN;
    dpx=2/(nx-1);
    parfor i=1:nx                      %Setting up the x momentum array list
        Px(i)=(i-nx/2-1/2)*dpx;
    end
    ny=200;
    Py(ny)=NaN;
    dpy=2/(ny-1);
    parfor i=1:ny                      %Do the same for the y momentum
        Py(i)=(i-ny/2-1/2)*dpy;    
    end

    %Momentum spectra
    Mx(nx,ny)=NaN;
    My(nx,ny)=NaN;
    P(nx,ny)=NaN;
    tic

    parfor ix=1:nx
        for iy=1:ny

            %K=sqrt(Px(ix)^2+Py(iy)^2);
            %Mx(ix,iy)=sqrt(128/3)*sqrt(K)*((K-I1)/(-K-I1))^(I1/k)*sqrt(1+coth(pi/K));
            %
            Mx(ix,iy)=(sqrt(-1)*2^(3.5)*(2*Ip)^(5/4)/pi)*Px(ix)/(Px(ix)^2+Py(iy)^2+2*Ip)^3; %Approximated field-free dipole transition element
            My(ix,iy)=(sqrt(-1)*2^(3.5)*(2*Ip)^(5/4)/pi)*Py(iy)/(Px(ix)^2+Py(iy)^2+2*Ip)^3;
            jkM=0;
            for it=1:nt
                phase=0;
                for it2=it:nt %phase = Integral (t to tf) (Px(x)-Ax(t))^2/2 + (Py(y)-Ay(t))^2/2
                    %phase=phase+dt*( (Px(ix)-Ax(it2))^2/2+(Py(iy)-Ay(it2))^2/2+Ip +0.1*(F(it)+Fd(it)) );
                    phase=phase+dt*( (Px(ix)+Ax(it2))^2/2+(Py(iy)+Ay(it2))^2/2+Ip );    %Quasiclassical Action (Where is speed of light?)
                end
                %jkM = Integral (1 to tf) i*(Ex(t)*M(x,y)+Ey(t)*M(x,y))*e^(-i*phi)dt
                jkM=jkM+sqrt(-1)*dt*Ex(it)*Mx(ix,iy)*exp(-sqrt(-1)*phase);
                jkM=jkM+sqrt(-1)*dt*Ey(it)*My(ix,iy)*exp(-sqrt(-1)*phase);
                P(ix,iy)=jkM;
            end
        end
        ix
    end
    pamp = log(abs(P').^2);
    fig = figure('visible','off');
    imagesc(Px,Py,pamp);
    set(gca, 'YDir','normal');
    colorbar;
    xlabel('Px (a.u.)');
    ylabel('Py (a.u.)');
    axis([-1 1 -1 1]);
    title('Photoelectron Momentum Distribution');
    str1 = {strcat('$$\tau = ', num2str(round(tau,1)), '\ au$$'), ...
            strcat('$$T = ', num2str(round(T,1)), '\ au$$'), ...
            strcat('$$\lambda = ', num2str(round(lambda,1)),'\ nm$$')};
    text(-0.95,0.80,str1,'Interpreter','latex','BackgroundColor','yellow');
    str2 = {strcat('$$\phi_1 = ', num2str(round(phi1/pi,2)), '\ \pi$$'), ...
            strcat('$$\phi_2 = ', num2str(round(phi2/pi,2)), '\ \pi$$')};
    text(-0.95,-0.85,str2,'Interpreter','latex','BackgroundColor','yellow');
    str3 = {strcat('$$n_c = ', num2str(nc),'$$')};
    text(0.5,0.90,str3,'Interpreter','latex','BackgroundColor','yellow');
    str4 = {strcat('$$Polarization = ', polarization,'$$')};
    text(0.05,-0.90,str4,'Interpreter','latex','BackgroundColor','yellow');
    pbaspect([1 1 1]);
    
    %Frige peak separation
    Prob = abs(P(100,1:100)');
    [arlen, arwid] = size(Prob);    %Get array length of 1D array
    maxima = [];                    %Instantiate maxima list
    for i=2:arlen-2                 %Iterate over indices with indices = index-1 and index+1 present
       if Prob(i) > Prob(i-1) && Prob(i) > Prob(i+1)
           maxima = [maxima; i];    %If value at index i is greater than adjacent values, add index to maxima list
       end
    end
    maxima_vals = Prob(maxima);     %Values at maxima
    max1 = maxima(maxima_vals==max(maxima_vals));     %Index of largest maxima
    max2 = maxima(maxima_vals==max(maxima_vals(maxima_vals<max(maxima_vals))));   %Index of 2nd largest maxima
    
    delta_p = abs(Px(max1) - Px(max2)); %p separation
    delta_E = abs(Px(max1)^2-Px(max2)^2)/2*27.211;  %E separation
    strp = {strcat('$$\Delta p_F = ',num2str(round(delta_p,2)),'\ au$$'), ...
        strcat('$$\Delta E_F = ',num2str(round(delta_E,2)),'\ eV$$')};
    hold on
    plot([Px(max1),Px(max2)],[0,0],'Color','r','LineWidth',2);  %Draws line between peaks
    text((Px(max1)+Px(max2))/2,0.2,strp,'Interpreter','latex','BackgroundColor','cyan');

    saveas(fig,strcat(dir,'Momentum_Distribution_logAmp_',zeroes,num2str(i_phi),'.png'))

    
    %Energy radial distribution
    PAmp = abs(P').^2;
    %EofR = zeros(int32(round(nx^2+ny^2)),n_tau);
    nr = round(sqrt((nx/2)^2+(ny/2)^2));
    r = 1:nr;
    EAmpofR(nr,n_phi) = NaN;
    PofR(nr) = NaN;
    parfor ir=1:nr
        for ix=1:nx
            for iy=1:ny
                if ir == round(sqrt((ix-nx/2)^2+(iy-ny/2)^2))
                    str = strcat(num2str(ir),',', num2str(ix),',', ...
                        num2str(iy),',', num2str(sqrt(Px(ix)^2+Py(iy)^2)),'\r\n');
                    PofR(ir) = sqrt(Px(ix)^2+Py(iy)^2);
                    EAmpofR(ir,i_phi) = EAmpofR(ir,i_phi) + PAmp(ix,iy);
                end
            end
        end
    end
    %save(strcat(dir,'data.mat'),'EAmpofR','PofR','phi_list','nr','n_phi','dir');
    
    
end
toc

%%
EofR(nr) = NaN;
for i_phi=1:n_phi
    for ir=1:nr
       EofR(ir) = PofR(ir)^2/2;
       EAmpofR(ir,i_phi) = EAmpofR(ir,i_phi)/max(max(EAmpofR(:,:)));
    end
end
EAmpofR(:,n_phi)=0;
phi_list = round(phi_list./pi, 4);
EAmpofR = EAmpofR;

fig_E_phi2 = figure('visible','on');
contourf(phi_list,EofR,EAmpofR,50,'edgecolor','none');
axis([min(phi_list) max(phi_list) 0 max(EofR)]);
colorbar;
%caxis([0 1]);
xlabel('$$\Delta\phi = \phi_2-\phi_1 = \phi_2-0=\phi_2\ (in\ multiples\ of\ pi)$$','Interpreter','latex');
ylabel('Energy (a.u.)');
title('Photoelectron Energy Distribution');
saveas(fig_E_phi2,strcat(dir,'E_vs_phi_contour.png'))

%%

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
