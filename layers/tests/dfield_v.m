function [ v dv ] = dfield_v(lay,gi,robs,rsrc)
% [ v dv ] = dfield_v(lay,gi,robs,rsrc)
%
% Used by dfield_dcim, calculates scalar potential of a point charge.
%

ki = lay.freq * sqrt(lay.eps2 .* mu0);  % Source/observation layer hardcoded

dr = robs-rsrc;
[ dx dy dz ] = uncat(1,dr);
[ xobs yobs zobs ] = uncat(1,robs);
[ xsrc ysrc zsrc ] = uncat(1,rsrc);

imgz0 = zobs-zsrc;
imgz1 = 2*lay.z3-(zobs+zsrc); % Source/observation layer hardcoded
imgz2 = (zobs+zsrc)-2*lay.z2;
imgz3 = 2*(lay.z3-lay.z2)+(zobs-zsrc);
imgz4 = 2*(lay.z3-lay.z2)-(zobs-zsrc);

[ Kf0 dKf0 ] = sumdcimg(gi.Kf0.a, gi.Kf0.b, dx, dy, imgz0, ki);
[ Kf1 dKf1 ] = sumdcimg(gi.Kf1.a, -gi.Kf1.b, dx, dy, -imgz1, ki); % ! Minus
[ Kf2 dKf2 ] = sumdcimg(gi.Kf2.a, gi.Kf2.b, dx, dy, imgz2, ki);
[ Kf3 dKf3 ] = sumdcimg(gi.Kf3.a, gi.Kf3.b, dx, dy, imgz3, ki);
[ Kf4 dKf4 ] = sumdcimg(gi.Kf4.a, -gi.Kf4.b, dx, dy, -imgz4, ki); % ! Minus

v = [ Kf0 + Kf1 + Kf2 + Kf3 + Kf4 ];
dv = [ dKf0 + dKf1 + dKf2 + dKf3 + dKf4 ];
