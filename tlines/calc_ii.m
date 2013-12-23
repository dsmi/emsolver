function i = calc_ii(tl,zobs,iobs,zsrc,jsrc)
% i = calc_ii(tl,zobs,iobs,zsrc,jsrc)
%
% Part of the tlines calculator (see calc_tlines), given the unit
% shunt current source at an arbitrary point, calculate the current
% at an another point.
%   tl   - structure with transmission lines parameters and auxiliary
%          data as returned by calc_tlines.
%   zobs - observation coordinate.
%   iobs - observation tline index.
%   zsrc - source coordinate.
%   jsrc - source tline index.

if iobs<jsrc,
	% Find voltage at the left terminal of the source tline
	vt = calc_vterm_i(tl,zsrc,jsrc);
	% Next, find voltage at the right terminal of the observation tline
	vt = vt.*prod(tl.Tls(:,iobs+1:jsrc-1),2);
	% Finally, find current at the observation point
	i = vt.*calc_i_vterm(tl,zobs,iobs);
else if iobs>jsrc,
	i = -calc_vv(tl,zsrc,jsrc,zobs,iobs);
else
	n = iobs;
	ex1 = exp(-tl.k(:,n).*(2*tl.z(:,n+1)-(zobs+zsrc)));
	i1 = -tl.Ggr(:,n).*ex1;
	ex2 = exp(-tl.k(:,n).*((zobs+zsrc)-2*tl.z(:,n)));
	i2 = tl.Gls(:,n).*ex2;
	ex3 = exp(-tl.k(:,n).*(2*tl.d(:,n)+(zobs-zsrc)));
	i3 = tl.Gls(:,n).*tl.Ggr(:,n).*ex3;
	ex4 = exp(-tl.k(:,n).*(2*tl.d(:,n)-(zobs-zsrc)));
	i4 = -tl.Gls(:,n).*tl.Ggr(:,n).*ex4;
	ex0 = exp(-tl.k(:,iobs).*abs(zobs-zsrc));
	if zobs>zsrc,
		i0 = ex0;
	else
		i0 = -ex0;
	end
	i = i0 + (i1+i2+i3+i4)./(1-tl.Gls(:,n).*tl.Ggr(:,n).*tl.t(:,n));
	i = i/2;
end
end
