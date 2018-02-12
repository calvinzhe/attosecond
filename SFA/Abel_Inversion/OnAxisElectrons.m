%%%OnAxisElectrons.m
%Author: Davy Foote
%Version: 1.0
%Date: July 31, 2017

%This program takes a series of 130 Abel-inverted electron images from the
%experiment on July 14, 2017, and finds the on-axis electron distributions
%for each image.

%The goal in the future will be to modify the program to look at electron
%distributions like this as a function of angle. This will be explained in
%the comments throughout this program.

clear, clc
close all

C=open('AbelInvertedImages.mat');
C=C.C;

%%%Plot an example electron image. This image shows the electron population
%%%P(r,z), where the laser propagation direction is z and r^2=x^2+y^2. The
%%%inversion algorithm assumes cylindrical symmetry, i.e, no dependence on
%%%the azimuthal angle phi.

ExampleImage=squeeze(C(:,:,1));

figure, imagesc(ExampleImage), axis square, colormap gray

xlabel('r (pixels)');
ylabel('z (pixels)');

UR = flipud(ExampleImage);
LR = ExampleImage;
LL = fliplr(ExampleImage);
UL = flipud(LL);

TotalImage(500,500)=NaN;
TotalImage(1:250,1:250)=UL;
TotalImage(1:250,251:500)=UR;
TotalImage(251:500,1:250)=LL;
TotalImage(251:500,251:500)=LR;

figure, imagesc(TotalImage), axis square, colormap gray
xlabel('r (pixels)');
ylabel('z (pixels)');

%%%NOTE: Indices for a 2D matrix in Matlab are ordered (row#, column#).
%%%This means that in this figure, the first index is the vertical axis r 
%%%and the second index is the horizontal axis z. If you look at a point on
%%%the figure, it is labeled (X,Y), which means these indices are reversed 
%%%from the array indices.

%%%Define rho^2=r^2+z^2.
%%%Define the angle alpha between the xy plane and z axis such that 
%%%r=rho*cos(alpha) and z=rho*sin(alpha).
%%%Using the standard spherical coordinates polar angle theta,
%%%alpha=pi/2-theta.

P=zeros(250,130); %initialize the array
tau=zeros(1,130); %tau is the delay between peaks of the TPP

for i=1:130 %Step through each image in the file
    P(:,i)=C(:,1,i); 
    %%%All rows of the first column of C. This gives us P(rho) for z=0,
    %%%equivalently alpha=0.
    
    tau(i)=130+0.057*i;
    %%%Each image corresponds to a change in tau of 0.057 fs.
end


%Plot the results
figure, imagesc(tau, 1:250, P);
set(gca,'YDir','normal') %This sets the (0,0) pixel as the lower left hand corner
xlabel('\tau (fs)');
ylabel('\rho (px)');

title('P(\rho), \alpha=\pi/2');

%%%The hope would be to get a similar image for different values of alpha
%%%between 0 and \pi/2 rad. Also, we will have to convert rho to energy, 
%%%U \propto rho^2 in the ideal case, but this can be done later.

%{
for i = 1:130
    for x = 1:250
        for y = 1:250
            PofR(int32(round(sqrt(x^2+y^2))),i) = PofR(int32(round(sqrt(x^2+y^2))),i)+C(x,y,i);
        end
    end
end
%}

