clear
N = 10000;
t = linspace(-500,500,N);

A=1;
B=1;
phi1=0;
phi2=0;
tau = 200;
f=3*10^8/(800*10^7.5);
sigma_t = 100;

N_tau = 500;
tau(N_tau) = NaN;
for i_tau=1:N_tau
    tau(i_tau) = i_tau*700/N_tau;
end

G1(N,N_tau)=NaN;
G2(N,N_tau)=NaN;
F1(N,N_tau)=NaN;
F2(N,N_tau)=NaN;
E(N,N_tau)=NaN;
I(N,N_tau)=NaN;

for i_tau=1:N_tau
    for i=1:N
        G1(i,i_tau) = A*exp(-(t(i)-tau(i_tau)/2)^2/(2*sigma_t^2));
        G2(i,i_tau) = B*exp(-(t(i)+tau(i_tau)/2)^2/(2*sigma_t^2));
        F1(i,i_tau) = cos(2*pi*f*(t(i)-tau(i_tau)/2)+phi1)*G1(i,i_tau);
        F2(i,i_tau) = cos(2*pi*f*(t(i)+tau(i_tau)/2)+phi2)*G2(i,i_tau);
        E(i,i_tau) = F1(i,i_tau)+F2(i,i_tau);
        I(i,i_tau) = E(i,i_tau)^2;
    end
end
%%
parfor i_tau=1:N_tau
    fig = figure();
    set(gcf,'Visible','off');
    subplot(2,1,1);
    plot(t,I(:,i_tau))
    ylim([0,4]);
    xlim([-500,500]);
    xlabel('Time (fs)');
    title('Intensity');
    
    subplot(2,1,2)
    plot(t,G1(:,i_tau),'r')
    hold on
    plot(t,G2(:,i_tau),'b')
    hold on
    plot(t,-G1(:,i_tau),'r')
    hold on
    plot(t,-G2(:,i_tau),'b')
    hold on
    plot(t,F1(:,i_tau),'r')
    hold on
    plot(t,F2(:,i_tau),'b')
    xlabel('Time (fs)');
    ylim([-2,2]);
    xlim([-500,500]);
    title('Electric Field');
    
    tau_str=strcat('$$\tau = ',num2str(tau(i_tau)),'\ $$fs');
    annotation('textbox',[0.15,0.6,0.3,0.3],'String',tau_str,...
        'Interpreter','latex','FitBoxToText','on');
    
    if i_tau < 10
        zeros_char = '00';
    elseif i_tau < 100
        zeros_char = '0';
    else
        zeros_char = '';
    end
            
    str = strcat('./sim_figures/interference_tau_',zeros_char,num2str(i_tau),'.png');
    saveas(fig,str);
end






