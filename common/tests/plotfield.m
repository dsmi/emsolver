% This script computes field of a Hertzian dipole, plots it and, optionally,
% saves it to a data file so it can be plot on Data Explorer using 
% vis/dipole.net visual program.

addpath(genpath([ pwd, '/..' ]));

% Free space
eps = eps0;
mu = mu0;

% Angular frequency
freq = 2e10;

% Wavenumber
k = freq * sqrt(eps * mu);
wavelen = 2*pi/k

% Observation points
xgrid = linspace( -0.1, 0.1, 200);
ygrid = linspace( -0.1, 0.1, 200);
[ x, y ] = meshgrid(xgrid, ygrid);
z = 0.009*ones(size(x));
 
% Y-oriented (vertical) dipole.
dy = ones(size(z));
dx = dz = zeros(size(x));

% Length of the dipole
dl = 0.03*ones(size(x));

% Dipole field.
[ ex, ey, ez, hx, hy, hz ] = dfield(freq, eps, mu, -x, -y, -z, dx, dy, dz, dl);

% Now the image dipole - under the z=0 plane, which we assume to
% be an ideal ground plane.
%x1 = x2 = y1 = y2 = zeros(size(x));
%z1 = -dpz*ones(size(x)) + dplen/2;
%z2 = -dpz*ones(size(x)) - dplen/2;

% Field of the image dipole
%[ iex, iey, iez ] = dpfield(freq,eps,mu,x,y,z,x1,y1,z1,x2,y2,z2);
%ex = ex + iex;
%ey = ey + iey;
%ez = ez + iez;

% Save field to file to plot it in Data Explorer
%px = x(:).';
%pz = z(:).';
%rex = (real(ex))(:).';
%rez = (real(ez))(:).';
%save("data.txt", "px", "pz", "rex", "rez", "-ascii")

val = abs(real(hz));

%contour(x,y,val);
imagesc(x,y,val);
%quiver(x,y,real(ex),real(ey))
