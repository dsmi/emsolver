function [ efj, efm ] = dfield_dcim(lay,gi,robs,rsrc,l)
% [ efj, efm ] = dfield_dcim(lay,gi,robs,rsrc,l)
%
% Calculates field of a current dipole in layered media via the Greens
% function in mixed potential form. First, the Greens functions for
% the scalar and vector potentials (together with the correction term)
% are calculated, and then the fields are computed from potentials.
% The Greens functios are available in closed form in spectral domain only,
% and the spatial domain values are approximated by complex images.
%  Inputs:
%    lay  - structure describing the layered media parameters, see
%           lay_defaults.m
%    gi   - complex images calculated by mkimages
%    robs - observation point, 1-by-3
%    rsrc - source point, 1-by-3
%    l    - source dipole, 1-by-3
%  Output:
%    efj  - electric field due to an electric current dipole, 1-by-3
%    efm  - electric field due to a magnetic current dipole, 1-by-3
%     

% Vector potential
[ Gvv Gzx Gzy Gzz ] = dfield_a(lay,gi,robs,rsrc);

% Correction term
[ Cxz Cyz Czz ] = dfield_c(lay,gi,robs,rsrc);

% Scalar potentials of the charges
[ v1 dv1 ] = dfield_v(lay,gi,robs,rsrc-l./2);
[ dvx1 dvy1 dvz1 ] = uncat(ndims(dv1), shiftdim(dv1, 1));
[ v2 dv2 ] = dfield_v(lay,gi,robs,rsrc+l./2);
v2 = -v2;  % Charge at the end of the dipole is negative
dv2 = -dv2;
[ dvx2 dvy2 dvz2 ] = uncat(ndims(dv2), shiftdim(dv2, 1));

[ lx ly lz ] = uncat(1,l);

% Compute the electric field from the vector and scalar potentials
% and the correction term.
efj = zeros(3,1);
efj(1) =  -Gvv*lx                - ( dvx1 + dvx2 ) - Cxz*lz;
efj(2) =  -Gvv*ly                - ( dvy1 + dvy2 ) - Cyz*lz;
efj(3) = -(Gzx*lx+Gzy*ly+Gzz*lz) - ( dvz1 + dvz2 ) - Czz*lz;

% Derivative of the electric vector potential
[ dFvv dFzx dFzy dFzz ] = dfield_df(lay,gi,robs,rsrc);

dflx = dFvv*lx;
dfly = dFvv*ly;
dflz = dFzx*lx+dFzy*ly+dFzz*lz;

% Compute the electric field as a curl of the electric vector potential
epsi = lay.eps2;
efm = zeros(3,1);
efm(1) = -1/(j*lay.freq*epsi)*(dflz(2)-dfly(3));
efm(2) = -1/(j*lay.freq*epsi)*(dflx(3)-dflz(1));
efm(3) = -1/(j*lay.freq*epsi)*(dfly(1)-dflx(2));
