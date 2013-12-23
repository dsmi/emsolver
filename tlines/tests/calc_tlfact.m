function [ iv ] = calc_tlfact(tl, f0, f1, f2, f3, f4, iobs, zobs, zsrc)
% Utility function used in the tlines calculator tests, given the 
% transmission line Greens function in factorized form (see fact_vi,
% fact_vv, fact_iv, fact_ii) and source and observation coordinates
% computes the voltage/current value.
%   f0-f4 - Factors, see expression below
%   iobs  - source (and observation) tline index
%   zobs  - observation coordinate, just z in the expression below.
%   zsrc  - source coordinate, z' in the expression below.
% The result is calculated using the following formula:
%   iv = f0*exp(-k*abs(z-z')) + f1*exp(-k*(2*z(i+1)-(z+z'))) 
%       + f2*exp(-k*((z+z')-2*z(i))) + f3*exp(-k*(2*d(i)+(z-z')))
%       + f4*exp(-k*(2*d(i)-(z-z')))
%

k = tl.k(:,iobs);
d = tl.d(:,iobs);
zn = tl.z(:,iobs);
zn1 = tl.z(:,iobs+1);

iv = f0*exp(-k*abs(zobs-zsrc)) + f1*exp(-k*(2*zn1-(zobs+zsrc))) ...
  + f2*exp(-k*((zobs+zsrc)-2*zn)) + f3*exp(-k*(2*d+(zobs-zsrc))) ...
	   + f4*exp(-k*(2*d-(zobs-zsrc)));
