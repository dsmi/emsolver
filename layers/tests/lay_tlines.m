function [ tle tlh ] = lay_tlines(lay,kr,kz)
% [ tle tlh ] = lay_tlines(lay,kr,kz)
%
% Sets up chain of transmission lines, which are used when computing
% spectral-domain Green's function for the layered media.
% There are two chains being set up, 'e' which corresponds to the TM
% fields and 'h' which corresponds to TE fields.
% Both tle and tlh are the outputs of calc_tlines function called with
% the corresponding transmission line parameters.
%

freq = lay.freq;

% Free-space wavenumber.
k0 = freq * sqrt(eps0*mu0);

% Flag indicating that there is an ideal ground above the top
% layer (layer N). Otherwise the first layer extends to infinity.
top_gnd = 0;

% Flag indicating that there is an ideal ground below the bottom
% layer (layer 1). Otherwise the bottom layer extends to infinity.
bottom_gnd = 0;

% Electric permittivity of the layers. Vector of length N.
lay_eps = [ lay.eps1 lay.eps2 lay.eps3 ];

% Magnetic permeability of the layers. Vector of length N.
lay_mu = [ mu0 mu0 mu0 ];

% Z-coordinate of the layers interfaces. Vector of length N+1, where N is
% the number of layers.
lay_z = [ lay.z2-1 lay.z2 lay.z3 lay.z3+1 ];

% Perfect electric conductor at the bottom?
if isinf(lay.eps1),
	bottom_gnd = 1;
	lay_eps(1) = [];
	lay_mu(1) = [];
	lay_z(1) = [];
	kz(:,1) = [];
end

% Perfect electric conductor at the top?
if isinf(lay.eps3),
	top_gnd = 1;
	lay_eps(end) = [];
	lay_mu(end) = [];
	lay_z(end) = [];
	kz(:,end) = [];
end

% Number of kr points
npt = length(kr);

Gls1_e = GgrN_e = Gls1_h = GgrN_h = 0;

% PEC : Ze = 0; Zh = 0, corresponds to short termination.
if bottom_gnd,
	Gls1_e = Gls1_h = -1;
end

if top_gnd,
	GgrN_e = GgrN_h = -1;
end

% This is common for both TM and TE tlines
tl_z = repmat(lay_z, npt, 1);
tl_k = kz*j; % Notice j

% Setup tlines associated with TM fields ('e' tlines)
Z0 = kz./repmat(freq*lay_eps, npt, 1);
Gls1 = repmat(Gls1_e, npt, 1);
GgrN = repmat(GgrN_e, npt, 1);

tle = calc_tlines(tl_z, Z0, tl_k, Gls1, GgrN);

% Setup tlines associated with TE fields ('h' tlines)
Z0 = repmat(freq*lay_mu, npt, 1)./kz;
Gls1 = repmat(Gls1_h, npt, 1);
GgrN = repmat(GgrN_h, npt, 1);

tlh = calc_tlines(tl_z, Z0, tl_k, Gls1, GgrN);
