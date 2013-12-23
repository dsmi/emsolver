function [ efj efm ] = dfield_dyg(lay,robs,iobs,rsrc,jsrc,l)
% [ efj efm ] = dfield_dyg(lay,robs,iobs,rsrc,jsrc,l)
%
% Calculates field of a current dipole in layered media via the dyadic
% Greens function. First, the function builds the dyadic Greens function
% for the given source and observation positions. It is available in closed
% form in spectral domain only, so it is first calulated in spectral domain
% with the help of the transmission lines calculator, and then transformed
% to spatial domain by evaluating the Sommerfeld integrals numerically.
% After the Greens function is available, the fields are calculated as:
%    E=<Gej;J>+<Gem;M>
%    H=<Ghj;J>+<Ghm;M>
% Currently it calculates the electric field due to electric and magnetic
% currents only.
% Works fairly slow, used in tests only. Does not work properly if z=z' or
% rho=rho'. 
%  Inputs:
%    lay  - structure describing the layered media parameters, see
%           lay_defaults.m
%    robs - observation point, 1-by-3
%    iobs - observation layer
%    rsrc - source point, 1-by-3
%    jsrc - source layer
%    l    - source dipole, 1-by-3
%  Output:
%    efj   - electric field due to an electric current dipole, 1-by-3
%    efm   - electric field due to a magnetic current dipole, 1-by-3
%

freq = lay.freq;

[ xobs yobs zobs ] = uncat(1,robs);
[ xsrc ysrc zsrc ] = uncat(1,rsrc);

% Angle between rho and the x-axis
phi = atan2(ysrc-yobs, xsrc-xobs);

% Length of projection of the observation-to-source vector
% to the X-Y plane
rho = sqrt((xsrc-xobs).*(xsrc-xobs)+(ysrc-yobs).*(ysrc-yobs));

% Permittivities of the source and observation layers.
lay_eps = [ lay.eps1 lay.eps2 lay.eps3 ];
epsi = lay_eps(iobs);
epsj = lay_eps(jsrc);
mui = mu0;
muj = mu0;

% If the first layer is perfectly conducting, then lay_tlines omits
% the first tline
if isinf(lay.eps1),
	iobs = iobs - 1;
	jsrc = jsrc - 1;
end

vih = @(kr,kz) calc_vi(lay_tlines_h(lay,kr,kz),zobs,iobs,zsrc,jsrc);

vih_minus_vie = @(kr,kz) ...
    calc_vi(lay_tlines_h(lay,kr,kz),zobs,iobs,zsrc,jsrc) ...
  - calc_vi(lay_tlines_e(lay,kr,kz),zobs,iobs,zsrc,jsrc);

vih_kr = @(kr,kz) vih(kr,kz).*kr;
vih_minus_vie_kr = @(kr,kz) vih_minus_vie(kr,kz).*kr;

s0_vih = eval_sommerf(lay, 0, vih_kr, rho);
s0_vih_minus_vie_kr2 = eval_sommerf(lay, 0, vih_minus_vie_kr, rho);
s2_vih_minus_vie = eval_sommerf(lay, 2, vih_minus_vie_kr, rho);

% In the spectral domain:
%  Gxx = vih+j*kx*j*kx*(vih-vie)/(kr*kr)
%  Gxx = vih+j*ky*j*ky*(vih-vie)/(kr*kr)
%  Gxy = i*kx*i*ky*(vih-vie)/kr*kr
%  Gyx = i*kx*i*ky*(vih-vie)/kr*kr
%
% Spectral->spatial domain:
%            Finv{G} = S0{F}
%      Finv{-j*kx*G} = -cos(phi)*S1{G}
%      Finv{-j*ky*G} = -sin(phi)*S1{G}
%  Finv{j*kx*j*kx*G} = 1/2*(cos(2*phi)*S2{G}-S0{kr*kr*F})
%  Finv{j*ky*j*ky*G} = -1/2*(cos(2*phi)*S2{G}+S0{kr*kr*F})
%  Finv{j*kx*j*ky*G} = 1/2*sin(2*phi)*S2{G}
%
Gxx = -(s0_vih+1/2*(cos(2*phi)*s2_vih_minus_vie-s0_vih_minus_vie_kr2));
Gyy = -(s0_vih-1/2*(cos(2*phi)*s2_vih_minus_vie+s0_vih_minus_vie_kr2));
Gxy = Gyx = -1/2*sin(2*phi)*s2_vih_minus_vie;

iie = @(kr,kz) calc_ii(lay_tlines_e(lay,kr,kz),zobs,iobs,zsrc,jsrc);
iie_kr2 = @(kr,kz) iie(kr,kz).*kr.*kr;

s1_iie = eval_sommerf(lay, 1, iie_kr2, rho);
Gzx = -cos(phi)*s1_iie/(j*freq*epsi);
Gzy = -sin(phi)*s1_iie/(j*freq*epsi);

ive = @(kr,kz) calc_iv(lay_tlines_e(lay,kr,kz),zobs,iobs,zsrc,jsrc);
ive_kr3 = @(kr,kz) ive(kr,kz).*kr.*kr.*kr;

s0_ive_kr2 = eval_sommerf(lay, 0, ive_kr3, rho);
Gzz = 1/(j*freq*epsj)*1/(j*freq*epsi)*s0_ive_kr2;

vve = @(kr,kz) calc_vv(lay_tlines_e(lay,kr,kz),zobs,iobs,zsrc,jsrc);
vve_kr2 = @(kr,kz) vve(kr,kz).*kr.*kr;

s1_vve = eval_sommerf(lay, 1, vve_kr2, rho);
Gxz = -cos(phi)*s1_vve/(j*freq*epsj);
Gyz = -sin(phi)*s1_vve/(j*freq*epsj);

[ lx ly lz ] = uncat(1,l);

efj = zeros(3,1);
efj(1) = Gxx*lx+Gxy*ly+Gxz*lz;
efj(2) = Gyx*lx+Gyy*ly+Gyz*lz;
efj(3) = Gzx*lx+Gzy*ly+Gzz*lz;

% Electric field due to the magnetic current in the spectral domain:
%   Fxx = -j*kx*j*ky/kr^2*(vve-vvh)
%   Fyy = -j*kx*j*ky/kr^2*(vvh-vve)
%   Fxy = j*kx*j*kx/kr^2*vve+j*ky*j*ky/kr^2*vvh
%   Fyx = -j*ky*j*ky/kr^2*vve-j*kx*j*kx/kr^2*vvh
%   Fzx = -ky/(freq*eps)*ive
%   Fzy = kx/(freq*eps)*ive
%   Fxz = ky/(freq*mu')*vih
%   Fyz = -kx/(freq*mu')*vih

vvh_minus_vve = @(kr,kz) ...
    calc_vv(lay_tlines_h(lay,kr,kz),zobs,iobs,zsrc,jsrc) ...
  - calc_vv(lay_tlines_e(lay,kr,kz),zobs,iobs,zsrc,jsrc);

vvh_minus_vve_kr = @(kr,kz) vvh_minus_vve(kr,kz).*kr;

s2_vvh_minus_vve = eval_sommerf(lay, 2, vvh_minus_vve_kr, rho);

Fxx = 1/2*sin(2*phi)*s2_vvh_minus_vve;
Fyy = -1/2*sin(2*phi)*s2_vvh_minus_vve;

vvh = @(kr,kz) calc_vv(lay_tlines_h(lay,kr,kz),zobs,iobs,zsrc,jsrc);

vvh_minus_vve = @(kr,kz) ...
    calc_vv(lay_tlines_h(lay,kr,kz),zobs,iobs,zsrc,jsrc) ...
  - calc_vv(lay_tlines_e(lay,kr,kz),zobs,iobs,zsrc,jsrc);

vvh_kr = @(kr,kz) vvh(kr,kz).*kr;
vvh_minus_vve_kr = @(kr,kz) vvh_minus_vve(kr,kz).*kr;

s0_vvh = eval_sommerf(lay, 0, vvh_kr, rho);
s0_vvh_minus_vve_kr2 = eval_sommerf(lay, 0, vvh_minus_vve_kr, rho);
s2_vvh_minus_vve = eval_sommerf(lay, 2, vvh_minus_vve_kr, rho);

Fyx = s0_vvh-1/2*(cos(2*phi)*s2_vvh_minus_vve+s0_vvh_minus_vve_kr2);
Fxy = -s0_vvh-1/2*(cos(2*phi)*s2_vvh_minus_vve-s0_vvh_minus_vve_kr2);

ive = @(kr,kz) calc_iv(lay_tlines_e(lay,kr,kz),zobs,iobs,zsrc,jsrc);
ive_kr2 = @(kr,kz) ive(kr,kz).*kr.*kr;

s1_ive = eval_sommerf(lay, 1, ive_kr2, rho);

Fzx = sin(phi)*s1_ive/(j*freq*epsi);
Fzy = -cos(phi)*s1_ive/(j*freq*epsi);

vih = @(kr,kz) calc_vi(lay_tlines_h(lay,kr,kz),zobs,iobs,zsrc,jsrc);
vih_kr2 = @(kr,kz) vih(kr,kz).*kr.*kr;

s1_vih = eval_sommerf(lay, 1, vih_kr2, rho);

Fxz = -sin(phi)*s1_vih/(j*freq*mu0);
Fyz = cos(phi)*s1_vih/(j*freq*mu0);

efm = zeros(3,1);
efm(1) = Fxx*lx+Fxy*ly+Fxz*lz;
efm(2) = Fyx*lx+Fyy*ly+Fyz*lz;
efm(3) = Fzx*lx+Fzy*ly;
