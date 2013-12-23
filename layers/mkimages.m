function gi = mkimages(lay)
% gi = mkimages(lay)
% 
% Approximate spatial domain Green's function with complex images.
%
% The input parameter is the structure with the layered media parameters,
% for description of the fields please see lay_default.m
%
% The results of the calculations are the potential Green's functions,
% which can be used to compute electric and magnetic fields:
%   E = -j*freq*mu0*<Ga;J> + 1/(j*freq*eps0)*grad(<Kf;div'(J)> + <Cf;J>)
%   H =
%
% Ga is a dyad with the following nonzero elements:
%
%  | Gvv     0     0 |
%  |   0   Gvv     0 |
%  | Gzu   Gzu   Gzz |
%
% Kf and Cf are the scalars.
%
% For each of the elements, z dependency is factored out, which gives
% the expansion of the following form:
%   A = A0*exp(-j*kz*abs(z-z')) + A1*exp(-j*kz*(2*z(i+1)-(z+z'))) 
%      + A2*exp(-j*kz*((z+z')-2*z(i))) + A3*exp(-j*kz*(2*d(i)+(z-z')))
%      + A4*exp(-j*kz*(2*d(i)-(z-z')))
% The layered Green's functions are derived in spectral domain, to get
% the spatial domain expressions the A0-A4 factros are approximated as
% a sum of exponentials, and transformed to the spatial domain using
% Sommerfeld identity - this is the discrete complex image method.
%
% Each component of the spectral domain Green's function is apprximated
% with a finite sum of exponentials of the following form.
%  G = sum(a.*exp(-b.*kz*j))/(j*kz)
%
% The fields of the output structure are:
%   Gvv0 Gvv1 Gvv2 Gvv3 Gvv4
%   Gzz0 Gzz1 Gzz2 Gzz3 Gzz4
%   Gzu1 Gzu2 Gzu3 Gzu4
%   Kf0 Kf1 Kf2 Kf3 Kf4
%   Cf1 Cf2 Cf3 Cf4
%   Fvv0 Fvv1 Fvv2 Fvv3 Fvv4
%   Fzz0 Fzz1 Fzz2 Fzz3 Fzz4
%   Fzu1 Fzu2 Fzu3 Fzu4
%

% Start from the quasi-static images
%gi = mkquasi(lay);

% instead of the quasistatic images just init all fields with the empty arrays
gs = gspec_calc(lay, 0, calc_kz(lay, 0)); % just to obtain the field names
gi = struct();
fldnames = fieldnames(gs);
for fnidx = 1:numel(fldnames)
	fn = fldnames{fnidx};
	ab.a = [];
	ab.b = [];
	gi = setfield(gi, fn, ab);
end

% Number of samples for both first and second paths
ns = 301;

k = lay.freq * sqrt(lay.eps2*mu0); % Source/observation layer hardcoded
jobs = 2;                          % Source/observation layer hardcoded

% Determine the integration path parameters
[ T01, T02 ] = calcppar(lay);

dt1 = T01/ns;
dt2 = T02/ns;

% First part of the path - straight
[ kr, kz, t ] = mkipath1(lay, T01, T02, ns);

% Compute the spectral domain values.
gs1 = gspec_calc(lay, kr, kz);

kzj = kz(:,jobs);

% Functions which transform the coefficients of the exponentials
fa1 = @(at, bt) at.*exp(-T02.*bt);
fb1 = @(at, bt) -bt/k;

% Approximate with complex images.
fldnames = fieldnames(gs1);  
for fnidx = 1:numel(fldnames)
	fn = fldnames{fnidx};
	specv = getfield(gs1, fn); % spectral domain values
	ab = getfield(gi, fn); % quasi-static images
	[ ab.a ab.b ] = appci(specv.*kzj*j, dt1, kzj, fa1, fb1, ab.a, ab.b);
	gi = setfield(gi, fn, ab);
end 

% Second part of the path - curved.
[ kr, kz, t ] = mkipath2(lay, T02, ns);

% Spectral domain values for the second part of the path
gs2 = gspec_calc(lay, kr, kz);

kzj = kz(:,jobs);

% Functions which transform the coefficients of the exponentials
fb2 = @(at, bt) bt*T02/(k*(j-T02));
fa2 = @(at, bt) at.*exp(k*j*fb2(at, bt));

% Complex images for the second part
fldnames = fieldnames(gs2);  
for fnidx = 1:numel(fldnames)
	fn = fldnames{fnidx};
	specv = getfield(gs2, fn); % spectral domain values
	ab = getfield(gi, fn); % images from the previous step
	[ ab.a ab.b ] = appci(specv.*kzj*j, dt2, kzj, fa2, fb2, ab.a, ab.b);
	gi = setfield(gi, fn, ab);
end 

% drop images with negligible contribution
%for fnidx = 1:numel(fldnames),
%	fn = fldnames{fnidx};
%	ab = getfield(gi, fn);
%	[ ab.a ab.b ] = dropimg(ab.a, ab.b);
%	gi = setfield(gi, fn, ab);
%end
