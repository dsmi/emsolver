function [ dFvv dFzx dFzy dFzz ] = dfield_df(lay,gi,robs,rsrc)
% [ dFvv dFzx dFzy dFzz ] = dfield_df(lay,gi,robs,rsrc)
%
% Used by dfield_dcim, calculates derivatives of the electric vector
% potential of a pointwise magnetic current.
%

ki = lay.freq * sqrt(lay.eps2 .* mu0);  % Source/observation layer hardcoded

dr = robs-rsrc;
[ dx dy dz ] = uncat(1,dr);
[ xobs yobs zobs ] = uncat(1,robs);
[ xsrc ysrc zsrc ] = uncat(1,rsrc);

imgz0 = zobs-zsrc;
imgz1 = 2*lay.z3-(zobs+zsrc);
imgz2 = (zobs+zsrc)-2*lay.z2;
imgz3 = 2*(lay.z3-lay.z2)+(zobs-zsrc);
imgz4 = 2*(lay.z3-lay.z2)-(zobs-zsrc);

[ Fvv0 dFvv0 ] = sumdcimg(gi.Fvv0.a, gi.Fvv0.b, dx, dy, imgz0, ki);
[ Fvv1 dFvv1 ] = sumdcimg(gi.Fvv1.a, -gi.Fvv1.b, dx, dy, -imgz1, ki); % ! Minus
[ Fvv2 dFvv2 ] = sumdcimg(gi.Fvv2.a, gi.Fvv2.b, dx, dy, imgz2, ki);
[ Fvv3 dFvv3 ] = sumdcimg(gi.Fvv3.a, gi.Fvv3.b, dx, dy, imgz3, ki);
[ Fvv4 dFvv4 ] = sumdcimg(gi.Fvv4.a, -gi.Fvv4.b, dx, dy, -imgz4, ki); % ! Minus

[ Fzz0 dFzz0 ] = sumdcimg(gi.Fzz0.a, gi.Fzz0.b, dx, dy, imgz0, ki);
[ Fzz1 dFzz1 ] = sumdcimg(gi.Fzz1.a, -gi.Fzz1.b, dx, dy, -imgz1, ki); % ! Minus
[ Fzz2 dFzz2 ] = sumdcimg(gi.Fzz2.a, gi.Fzz2.b, dx, dy, imgz2, ki);
[ Fzz3 dFzz3 ] = sumdcimg(gi.Fzz3.a, gi.Fzz3.b, dx, dy, imgz3, ki);
[ Fzz4 dFzz4 ] = sumdcimg(gi.Fzz4.a, -gi.Fzz4.b, dx, dy, -imgz4, ki); % ! Minus

[ Fzu1 dFzu1 d2Fzu1 ] = sumdcimg(gi.Fzu1.a, gi.Fzu1.b, dx, dy, imgz1, ki);
dFzx1 = -1/j*shiftdim(d2Fzu1(1,:,:), 1); % Notice minus - this is because
dFzy1 = -1/j*shiftdim(d2Fzu1(2,:,:), 1); % differentiation by dx/dy in spatial
                                     % domain corresponds to multiplication
                                     % by -j*kx / -j*ky in spectral domain.

% Notice minus - because we need derivative with respect to zobs
[ Fzu2 dFzu2 d2Fzu2 ] = sumdcimg(gi.Fzu2.a, -gi.Fzu2.b, dx, dy, -imgz2, ki);
dFzx2 = -1/j*shiftdim(d2Fzu2(1,:,:), 1);
dFzy2 = -1/j*shiftdim(d2Fzu2(2,:,:), 1);

[ Fzu3 dFzu3 d2Fzu3 ] = sumdcimg(gi.Fzu3.a, gi.Fzu3.b, dx, dy, imgz3, ki);
dFzx3 = -1/j*shiftdim(d2Fzu3(1,:,:), 1);
dFzy3 = -1/j*shiftdim(d2Fzu3(2,:,:), 1);

% Notice minus - because we need derivative with respect to zobs
[ Fzu4 dFzu4 d2Fzu4 ] = sumdcimg(gi.Fzu4.a, -gi.Fzu4.b, dx, dy, -imgz4, ki);
dFzx4 = -1/j*shiftdim(d2Fzu4(1,:,:), 1);
dFzy4 = -1/j*shiftdim(d2Fzu4(2,:,:), 1);

dFvv = dFvv0+dFvv1+dFvv2+dFvv3+dFvv4;
dFzz = dFzz0+dFzz1+dFzz2+dFzz3+dFzz4;
dFzx =       dFzx1+dFzx2+dFzx3+dFzx4;
dFzy =       dFzy1+dFzy2+dFzy3+dFzy4;
