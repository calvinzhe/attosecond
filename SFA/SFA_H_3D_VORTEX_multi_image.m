%SFA_H_3D_VORTEX_multi-image.m: SFA, real hydrogen atom, Feb20 2017
%works
clear,clc,tic
close all
I1=sqrt(-1);
opengl software;

%Two oppositely circularly polarized time-delayed as pulses: PRL 115 113004 2015
nt=500;
W=3/27.2;    %0.1103 angular frequency => Period = 2*pi/W = 56.9675 => F = 0.0176
lambda = (3*10^8)/(W/2*pi)*(2.419*10^-17)/10^-9;  %wavelength in nm
T=6*2*pi/W;  %Pulse width = N*Period = 6*[2*3.14/(3/27.2)] = 341.8 | e^(-T^2/T^2) = e^(-1) = 0.3679
%T=2000; %pulse width in atomic units
tau=T;       %Pulse separation
phi1 = 0;
phi2 = 0;
E0=sqrt(3.5e14/3.5e16);
Ip=12.13/27.2;  %He=0.9037 for He
ksi=1;
dt=8*T/(nt-1);   %8*341.8/(100-1) = 27.62 Time increment for loop iteration
t(nt)=NaN;
td(nt)=NaN;
Ex(nt)=NaN;
Ey(nt)=NaN;
F(nt)=NaN;
Fd(nt)=NaN;



for i_phi2=1:13
    phi2 = (i_phi2-1)*pi/12;
    i_phi2

    for i_tau=1:6
        i_plot = rem(i_tau, 3);
        if i_plot == 1
            fig = figure;
            set(fig, 'Position',[0 0 1280 420]);
            i_tau_3set = idivide(i_tau, int32(3))+1;
        elseif i_plot == 0
            i_plot = 3;
        end
        tau = (i_tau-1)*0.75*T;
        for i=1:nt
            t(i) =(i-(nt-1)/2)*dt+tau/2;     %Time variable for 1st pulse: [-29,70]dt increments of dt or [-800.98,1933.4] w/ center at 566.21
            td(i)=(i-(nt-1)/2)*dt-tau/2;     %Time variable for the 2nd (delayed) pulse: [-69,30]dt increments of dt or [-1905.8,816] w/ center at -544.9 | center separation 1111.11 = 3.2507*T
            F(i)=exp(-t(i)^2/T^2);      %Envelop function for 1st pulse, Gaussian
            Fd(i)=exp(-td(i)^2/T^2);    %Envelop fuction for the 2nd (delayed) pulse, Gaussian
            Ex(i)=   F(i)*E0    /(1+ksi^2)^0.5*cos(1.0*W*t(i)+phi1) + 1*Fd(i)*E0    /(1+ksi^2)^0.5*cos(-1.0*W*td(i)+phi2);
            Ey(i)=   F(i)*E0*ksi/(1+ksi^2)^0.5*sin(1.0*W*t(i)+phi1) + 1*Fd(i)*E0*ksi/(1+ksi^2)^0.5*sin(-1.0*W*td(i)+phi2);
        end
        plotrange = linspace(min(td),max(t),nt);
        
        %{
        subplot(5,6,i_plot);
        plot(plotrange,E0*F/sqrt(2),'b', ...
                plotrange,E0*Fd/sqrt(2),'r', ...
                plotrange,sqrt(Ex.^2+Ey.^2),'g');
        legend('1st Pulse Envelope','2nd Pulse Envelope','E-field Magnitude');
        xlabel('time');
        
        subplot(5,6,i_plot+1);
        plot3(plotrange,Ex,Ey);
        xlabel('time');
        ylabel('Ex');
        zlabel('Ey');
        set(gca, 'YDir','reverse');
        %}
        

        %Vector potential
        % E(t) = -dA(t)/dt  =>  A(t) = -Integral (t0 to t) E(t')dt'
        Ax(nt)=NaN;
        Ay(nt)=NaN;
        for i=1:nt
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
        nx=100;
        Px(nx)=NaN;
        dpx=3/(nx-1);
        for i=1:nx                      %Setting up the x momentum array list
            Px(i)=(i-nx/2-1/2)*dpx;     %Px in [-49.5,49.5] in increments of dpx=1/33
        end
        ny=100;
        Py(ny)=NaN;
        dpy=3/(ny-1);
        for i=1:ny                      %Do the same for the y momentum
            Py(i)=(i-ny/2-1/2)*dpy;    
        end

        %Momentum spectra
        Mx(nx,ny)=NaN;
        My(nx,ny)=NaN;
        P(nx,ny)=NaN;
        tic
        for ix=1:nx
            for iy=1:ny

                %K=sqrt(Px(ix)^2+Py(iy)^2);
                %Mx(ix,iy)=sqrt(128/3)*sqrt(K)*((K-I1)/(-K-I1))^(I1/k)*sqrt(1+coth(pi/K));
                %
                Mx(ix,iy)=(sqrt(-1)*2^(3.5)*(2*Ip)^(5/4)/pi)*Px(ix)/(Px(ix)^2+Py(iy)^2+2*Ip)^3;
                My(ix,iy)=(sqrt(-1)*2^(3.5)*(2*Ip)^(5/4)/pi)*Py(iy)/(Px(ix)^2+Py(iy)^2+2*Ip)^3;
                jkM=0;
                for it=1:nt
                    phase=0;
                    for it2=it:nt %phase = Integral (t to tf) (Px(x)-Ax(t))^2/2 + (Py(y)-Ay(t))^2/2
                        %phase=phase+dt*( (Px(ix)-Ax(it2))^2/2+(Py(iy)-Ay(it2))^2/2+Ip +0.1*(F(it)+Fd(it)) );
                        phase=phase+dt*( (Px(ix)+Ax(it2))^2/2+(Py(iy)+Ay(it2))^2/2+Ip );
                    end
                    %jkM = Integral (1 to tf) i*(Ex(t)*M(x,y)+Ey(t)*M(x,y))*e^(-i*phi)dt
                    jkM=jkM+sqrt(-1)*dt*Ex(it)*Mx(ix,iy)*exp(-sqrt(-1)*phase);
                    jkM=jkM+sqrt(-1)*dt*Ey(it)*My(ix,iy)*exp(-sqrt(-1)*phase);
                end
                P(ix,iy)=jkM;
            end
            ix
        end
        toc/60;
        
        subplot(1,3,i_plot);
        contour(Px,Py,abs(P').^2,50);
        xlabel('Px (a.u.)');
        ylabel('Py (a.u.)');
        axis([-1 1 -1 1]);
        title('Photoelectron Momentum Distribution');
        str1 = {strcat('$$\tau = ', num2str(round(tau,1)), '\ au$$'), ...
                strcat('$$T = ', num2str(round(T,1)), '\ au$$'), ...
                strcat('$$\lambda = ', num2str(round(lambda,1)),'\ nm$$')};
        text(-0.95,0.75,str1,'Interpreter','latex','BackgroundColor','yellow');
        str2 = {strcat('$$\phi_1 = ', num2str(round(phi1/pi,2)), '\ \pi$$'), ...
                strcat('$$\phi_2 = ', num2str(round(phi2/pi,2)), '\ \pi$$')};
        text(-0.95,-0.8,str2,'Interpreter','latex','BackgroundColor','yellow');
        pbaspect([1 1 1])
        
        if i_plot == 3
            saveas(fig,strcat('./CCP_tau_phi12_scan/p_vs_tau_-_phi2_',num2str(round(phi2/pi,2)),'_set_',num2str(i_tau_3set),'.png'));
        end
    end
end


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
