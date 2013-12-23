function test_sommerf
% Obtain the spatial domain values of the layered Green's functions
% in two ways: by evaluating the Sommerfeld integral numerically
% and by applying the Sommerfeld identity to the approximating
% exponentials; and compare the results.
%

% Angular frequency
lay.freq = 1e10;

% Setup the layers stackup
lay.z2   = 1e-3;
lay.z3   = 3e-3;
lay.eps1 = eps0*4;
lay.eps2 = eps0*10;
lay.eps3 = eps0*2;

% This call generates the complex images which approximate the
% spatial domain layered Green's function.
gi = mkimages(lay);

ki = lay.freq * sqrt(lay.eps2 .* mu0);  % Source/observation layer hardcoded

% Source and observation positions used in the tests
zsrc = 2.5e-3;  % z', source height
zobs = 1.5e-3;  % z, observation height
xsrc = 2.5e-3;  % x', source position
xobs = -7.0e-3; % x, observation position
ysrc = -3.2e-5; % y', source position
yobs = 4.0e-5;  % y, observation position

% Notice! When calculating dx and dy, we subtract the source coordinate
% from the observation coordinate, this is important because when
% computing Gzx and Gzy we need to calculate the derivatives by x and y
% unprimed, i.e. the observation point coordinates. Otherwise, the
% derivative changes sign.
dx = xobs-xsrc;
dy = yobs-ysrc;
dz = zobs-zsrc;

% Angle between rho and the x-axis
phi = atan2(ysrc-yobs, xsrc-xobs);

% Length of projection of the observation-to-source vector
% to the X-Y plane
rho = sqrt((xsrc-xobs).*(xsrc-xobs)+(ysrc-yobs).*(ysrc-yobs));

reltol = 2e-2; % Relative tolerance.

fldnames = fieldnames(gi);
for fnidx = 1:numel(fldnames)
	fldname = fldnames{fnidx};
	fg = @(kr, kz) getfield(gspec_calc(lay, kr, kz), fldname).*kr;
	G_sommerf = eval_sommerf(lay, 0, fg, rho);
	ab = getfield(gi, fldname); % images
	G_dcim = sumdcimg(ab.a, ab.b, dx, dy, 0, ki);
	assertEquals(G_sommerf, G_dcim, abs(G_sommerf)*reltol);
end 

% Next, test the approximation for the particular value of z and z'.
fex = @(kz, imgz) exp(-j*kz(:,2)*imgz);  % Source/observation layer hardcoded
imgz0 = abs(zobs-zsrc);
imgz1 = 2*lay.z3-(zobs+zsrc);
imgz2 = (zobs+zsrc)-2*lay.z2;
imgz3 = 2*(lay.z3-lay.z2)+(zobs-zsrc);
imgz4 = 2*(lay.z3-lay.z2)-(zobs-zsrc);

% Image positions for each of the fields
[ imz.Gvv0, imz.Gzz0, imz.Kf0 ] = deal(imgz0, imgz0, imgz0);
[ imz.Fvv0, imz.Fzz0 ] = deal(imgz0, imgz0);
[ imz.Gvv1, imz.Gvv2, imz.Gvv3, imz.Gvv4 ] = deal(imgz1, imgz2, imgz3, imgz4);
[ imz.Gzz1, imz.Gzz2, imz.Gzz3, imz.Gzz4 ] = deal(imgz1, imgz2, imgz3, imgz4);
[ imz.Gzu1, imz.Gzu2, imz.Gzu3, imz.Gzu4 ] = deal(imgz1, imgz2, imgz3, imgz4);
[ imz.Kf1,  imz.Kf2,  imz.Kf3,  imz.Kf4  ] = deal(imgz1, imgz2, imgz3, imgz4);
[ imz.Cf1,  imz.Cf2,  imz.Cf3,  imz.Cf4  ] = deal(imgz1, imgz2, imgz3, imgz4);
[ imz.Gzz1, imz.Gzz2, imz.Gzz3, imz.Gzz4 ] = deal(imgz1, imgz2, imgz3, imgz4);
[ imz.Fvv1, imz.Fvv2, imz.Fvv3, imz.Fvv4 ] = deal(imgz1, imgz2, imgz3, imgz4);
[ imz.Fzz1, imz.Fzz2, imz.Fzz3, imz.Fzz4 ] = deal(imgz1, imgz2, imgz3, imgz4);
[ imz.Fzu1, imz.Fzu2, imz.Fzu3, imz.Fzu4 ] = deal(imgz1, imgz2, imgz3, imgz4);

fldnames = fieldnames(gi);
for fnidx = 1:numel(fldnames)
	fldname = fldnames{fnidx};
	imgz = getfield(imz, fldname);
	fg = @(kr, kz) getfield(gspec_calc(lay, kr, kz), fldname).*fex(kz, imgz).*kr;
	G_sommerf = eval_sommerf(lay, 0, fg, rho);
	ab = getfield(gi, fldname); % images
	G_dcim = sumdcimg(ab.a, ab.b, dx, dy, imgz, ki);
	assertEquals(G_sommerf, G_dcim, abs(G_sommerf)*reltol);
end 

%
% Next, test approximation of the terms which depend on the relative
% position of source and observation point x-x' y-y'. The corresponding
% spectral domain expressions include linear dependence on kx and ky.
% For these terms, the spectral domain value with kx/ky dropped is
% approximated with exponentials, the exponentials get transformed
% to the spatial domain and then differentiated by x/y, because according
% to the Fourier transform properties multiplication by j*kx/j*ky in
% spectral domain corresponds to differentiation by x/y in spatial domain.
% 

fldnames = { 'Gzu1', 'Gzu2', 'Gzu3', 'Gzu4' };
for fnidx = 1:numel(fldnames)
	
	fldname = fldnames{fnidx};
	
	fg = @(kr, kz) getfield(gspec_calc(lay, kr, kz), fldname).*kr.*kr;
	Gzu_sommerf1 = eval_sommerf(lay, 1, fg, rho);
	Gzx_sommerf = -cos(phi)/(-j)*Gzu_sommerf1;
	Gzy_sommerf = -sin(phi)/(-j)*Gzu_sommerf1;

	ab = getfield(gi, fldname); % images
	[ Gzu_dcim dGzu_dcim  ] = sumdcimg(ab.a, ab.b, dx, dy, 0, ki);

	% Notice 1/j factor! - this factor appears because differentiation
	% in spatial domain corresponds to multiplication by j*kx/j*ky; while
	% we drop only kx/ky when approximating the corresponding term with
	% exponentials.
	Gzx_dcim = 1/j*shiftdim(dGzu_dcim(1,:), 1);
	Gzy_dcim = 1/j*shiftdim(dGzu_dcim(2,:), 1);
	assertEquals(Gzx_sommerf, Gzx_dcim, abs(Gzx_sommerf)*reltol);
	assertEquals(Gzy_sommerf, Gzy_dcim, abs(Gzy_sommerf)*reltol);	

end 
