% Plots the spectral domain Greens funtion

addpath(genpath([ pwd, '/../..' ]));

% Angular frequency
lay.freq = freq = 1e11;

% Setup the layers stackup
lay.z3   = 3e-3;
lay.z2   = 1e-3;
lay.eps3 = eps0*2;
lay.eps2 = eps0*10;
lay.eps1 = eps0-j*Inf;

regv = linspace(-1500, 1500, 50);
imgv = linspace(-500, 500, 50);
[ re im ] = meshgrid(regv, imgv);
kr = re + j*im;
kz = calc_kz(lay, kr);

kzj = kz(:,:,2); % Source/observation layer

gs = gspec_calc(lay, kr, kz);

%mesh(re, im, abs(gs.Cf2))

zsrc = 2.5e-3;  % z', source height
zobs = 1.5e-3;  % z, observation height

imgz0 = zobs-zsrc;
imgz1 = 2*lay.z3-(zobs+zsrc);
imgz2 = (zobs+zsrc)-2*lay.z2;
imgz3 = 2*(lay.z3-lay.z2)+(zobs-zsrc);
imgz4 = 2*(lay.z3-lay.z2)-(zobs-zsrc);

ex0 = exp(-j*kzj*imgz0);
ex1 = exp(-j*kzj*imgz1);
ex2 = exp(-j*kzj*imgz2);
ex3 = exp(-j*kzj*imgz3);
ex4 = exp(-j*kzj*imgz4);

Gvv = gs.Gvv1.*ex1 + gs.Gvv2.*ex2 + gs.Gvv3.*ex3 + gs.Gvv4.*ex4;

mesh(re, im, abs(gs.Gvv2))
