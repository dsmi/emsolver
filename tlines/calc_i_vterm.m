function i = calc_i_vterm(tl,zobs,iobs)
% i = calc_i_vterm(tl,zobs,iobs)
%
% Part of the tlines calculator (see calc_tlines), given unit voltage
% at the right terminal of a tline calculate current at an arbitrary
% point of this tline.
%   tl   - structure with transmission lines parameters and auxiliary
%          data as returned by calc_tlines.
%   zsrc - source coordinate.
%   jsrc - source tline index.

denom = 1 + tl.Gls(:,iobs).*tl.t(:,iobs);
ex1 = exp(-2*tl.k(:,iobs).*(zobs - tl.z(:,iobs)));
ex2 = exp(-tl.k(:,iobs).*(tl.z(:,iobs+1) - zobs));
i = -tl.Y0(:,iobs).*(1 - tl.Gls(:,iobs).*ex1).*ex2./denom;
