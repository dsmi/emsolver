function [ f0, f1, f2, f3, f4 ] = fact_ii(tl,iobs)
% [ f0, f1, f2, f3, f4 ] = fact_ii(tl,iobs)
%
% Part of the tlines calculator (see calc_tlines), calculates the
% factorized form of the transmission line Greens function which
% allows to find the current given the shunt current source.
% The function results are five multipliers, using these multipliers
% the current can be found according to the formula below.
%   i = f0*exp(-k*abs(z-z')) + f1*exp(-k*(2*z(i+1)-(z+z'))) 
%      + f2*exp(-k*((z+z')-2*z(i))) + f3*exp(-k*(2*d(i)+(z-z')))
%      + f4*exp(-k*(2*d(i)-(z-z')))
%  where:
%   i      - the resulting current
%   z      - observation coordinate
%   z'     - source coordinate coordinate
%   k      - transmission line propagation constant
%   z(i+1) - coordinate of the end of the transmission line
%   z(i)   - coordinate of the beginning of the transmission line
%   d(i)   - length of the transmission line; d(i) = z(i+1) - z(i)
%   f0     - multiplier which corresponds to the direct ray between
%            the souce and observation, returned by this function.
%            Notice that this Greens function has discontinuous f0
%            multiplier, it has different signs for the cases z>z'
%            and z<z'. The returned value corresponds to the z>z' case.
%   f1-f4 -  multipliers which correspond to the rays which undergo
%            partial reflections at the tline junctions, returned by
%            this function.
% Such factorized form is possible only if both source and observation
% points belong to the same tline.
%  Inputs:
%    tl   - structure with transmission lines parameters and auxiliary
%           data as returned by calc_tlines.
%    iobs - source/observation tline
%

n = iobs;

% Common multiplier
m14 = 1./((1-tl.Gls(:,n).*tl.Ggr(:,n).*tl.t(:,n))*2);

% Factors corresponding to the reflected rays
f1 = -tl.Ggr(:,n).*m14;
f2 = tl.Gls(:,n).*m14;
f3 = tl.Gls(:,n).*tl.Ggr(:,n).*m14;
f4 = -tl.Gls(:,n).*tl.Ggr(:,n).*m14;

% Factor corresponding to the direct ray
f0 = repmat(1/2, size(f1));
