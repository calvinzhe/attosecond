[r,phi] = meshgrid(linspace(0,10,100),linspace(0,2*pi,100));
X = r.*cos(phi);
Y = r.*sin(phi);

Eb = 1;     %Bound energy
tau = 1;    %Pulse separation
phi_12 = 0; %CEP difference
PHI = (r.^2/2 + Eb)*tau + phi_12;
tdp = cos(PHI/2-phi);
figure;
surf(X,Y,tdp,'EdgeColor','none');
view(2);