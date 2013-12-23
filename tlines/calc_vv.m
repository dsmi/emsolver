function v = calc_vv(tl,zobs,iobs,zsrc,jsrc)
% v = calc_vv(tl,zobs,iobs,zsrc,jsrc)
%
% Part of the tlines calculator (see calc_tlines), given the unit
% series voltage source at an arbitrary point, calculate the voltage
% at an another point.
%   tl   - structure with transmission lines parameters and auxiliary
%          data as returned by calc_tlines.
%   zobs - observation coordinate.
%   iobs - observation tline index.
%   zsrc - source coordinate.
%   jsrc - source tline index.

if iobs<jsrc,
	% Find voltage at the left terminal of the source tline
	vt = calc_vterm_v(tl,zsrc,jsrc);
	% Next, find voltage at the right terminal of the observation tline
	vt = vt.*prod(tl.Tls(:,iobs+1:jsrc-1),2);
	% Finally, find voltage at the observation point
	v = vt.*calc_v_vterm(tl,zobs,iobs);
else if iobs>jsrc,
	v = -calc_ii(tl,zsrc,jsrc,zobs,iobs);
else
	n = iobs;
	ex1 = exp(-tl.k(:,n).*(2*tl.z(:,n+1)-(zobs+zsrc)));
	v1 = tl.Ggr(:,n).*ex1;
	ex2 = exp(-tl.k(:,n).*((zobs+zsrc)-2*tl.z(:,n)));
	v2 = -tl.Gls(:,n).*ex2;
	ex3 = exp(-tl.k(:,n).*(2*tl.d(:,n)+(zobs-zsrc)));
	v3 = tl.Gls(:,n).*tl.Ggr(:,n).*ex3;
	ex4 = exp(-tl.k(:,n).*(2*tl.d(:,n)-(zobs-zsrc)));
	v4 = -tl.Gls(:,n).*tl.Ggr(:,n).*ex4;
	ex0 = exp(-tl.k(:,iobs).*abs(zobs-zsrc));
	if zobs>zsrc,
		v0 = ex0;
	else
		v0 = -ex0;
	end
	v = v0 + (v1+v2+v3+v4)./(1-tl.Gls(:,n).*tl.Ggr(:,n).*tl.t(:,n));
	v = v/2;
end
end
