function G = sum_sommerf(n, Gspec, kr, rho)
% Used by eval_sommerf - approximate Sommerfeld integral over the given
% range with sum.
% Inputs:
%   n     - order of the Sommerfeld integral
%   Gspec - spectral domain values multiplied by kr^(n+1) calculated
%           at a series of points along the integration path, vector
%   kr    - points where the Gspec has been computed, vector of the same
%           length as Gspec
%   rho   - Radial distance to the point where the spatial domain value
%           is to be evaluated.
% Output:
%   G     - the resulting value.
%

integrand = Gspec.*besselj(n,kr*rho);
dkr = kr(2:end) - kr(1:end-1);
G = 1/(2*pi)*sum((integrand(1:end-1)+integrand(2:end))*0.5.*dkr);
