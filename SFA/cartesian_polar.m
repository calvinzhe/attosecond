nx=200;
Px(nx)=NaN;
dpx=2/(nx-1);
parfor i=1:nx                      %Setting up the x momentum array list
    Px(i)=(i-nx/2-1/2)*dpx;     %Px in [-49.5,49.5] in increments of dpx=1/33
end
ny=200;
Py(ny)=NaN;
dpy=2/(ny-1);
parfor i=1:ny                      %Do the same for the y momentum
    Py(i)=(i-ny/2-1/2)*dpy;    
end

mat= zeros(nx,ny);
for ix=1:nx
    for iy=1:ny
        mat(ix,iy)=sqrt(Px(ix)^2+Py(iy)^2);
    end
end

nr = round(sqrt((nx/2)^2+(ny/2)^2));
mat2 = zeros(nr,1);
PofR = zeros(nr,1);
EofR = zeros(nr,1);
for ir=1:nr
    for ix=1:nx
        for iy=1:ny
            if ir == round(sqrt((ix-nx/2)^2+(iy-ny/2)^2))
                PofR(ir) = sqrt(Px(ix)^2+Py(iy)^2);
                EofR(ir) = PofR(ir)^2/2;
            end
        end
    end
end