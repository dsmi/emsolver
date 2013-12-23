function [ Cxz Cyz Czz ] = dfield_c(lay,gi,robs,rsrc)
% [ Cxz Cyz Czz ] = dfield_c(lay,gi,robs,rsrc)
%
% Used by dfield_dcim, calculates the correction terms associtated with
% a pointwise current.
%

ki = lay.freq * sqrt(lay.eps2 .* mu0);  % Source/observation layer hardcoded

dr = robs-rsrc;
[ dx dy dz ] = uncat(1,dr);
zobs = robs(3);
zsrc = rsrc(3);

imgz0 = zobs-zsrc;
imgz1 = 2*lay.z3-(zobs+zsrc);
imgz2 = (zobs+zsrc)-2*lay.z2;
imgz3 = 2*(lay.z3-lay.z2)+(zobs-zsrc);
imgz4 = 2*(lay.z3-lay.z2)-(zobs-zsrc);

[ Cf1 dCf1 ] = sumdcimg(gi.Cf1.a, -gi.Cf1.b, dx, dy, -imgz1, ki); % !Minus
[ Cf2 dCf2 ] = sumdcimg(gi.Cf2.a, gi.Cf2.b, dx, dy, imgz2, ki);
[ Cf3 dCf3 ] = sumdcimg(gi.Cf3.a, gi.Cf3.b, dx, dy, imgz3, ki);
[ Cf4 dCf4 ] = sumdcimg(gi.Cf4.a, -gi.Cf4.b, dx, dy, -imgz4, ki); % !Minus

dCf = dCf1 + dCf2 + dCf3 + dCf4;
[ Cxz Cyz Czz ] = uncat(1, dCf);
