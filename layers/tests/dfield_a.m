function [ Gvv Gzx Gzy Gzz  ] = dfield_a(lay,gi,robs,rsrc)
% [ Gvv Gzx Gzy Gzz  ] = dfield_a(lay,gi,robs,rsrc)
%
% Used by dfield_dcim, calculates vector potential of a pointwise current.
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

Gvv0 = sumdcimg(gi.Gvv0.a, gi.Gvv0.b, dx, dy, imgz0, ki);
Gvv1 = sumdcimg(gi.Gvv1.a, gi.Gvv1.b, dx, dy, imgz1, ki);
Gvv2 = sumdcimg(gi.Gvv2.a, gi.Gvv2.b, dx, dy, imgz2, ki);
Gvv3 = sumdcimg(gi.Gvv3.a, gi.Gvv3.b, dx, dy, imgz3, ki);
Gvv4 = sumdcimg(gi.Gvv4.a, gi.Gvv4.b, dx, dy, imgz4, ki);

Gzz0 = sumdcimg(gi.Gzz0.a, gi.Gzz0.b, dx, dy, imgz0, ki);
Gzz1 = sumdcimg(gi.Gzz1.a, gi.Gzz1.b, dx, dy, imgz1, ki);
Gzz2 = sumdcimg(gi.Gzz2.a, gi.Gzz2.b, dx, dy, imgz2, ki);
Gzz3 = sumdcimg(gi.Gzz3.a, gi.Gzz3.b, dx, dy, imgz3, ki);
Gzz4 = sumdcimg(gi.Gzz4.a, gi.Gzz4.b, dx, dy, imgz4, ki);

[ Gzu1 dGzu1  ] = sumdcimg(gi.Gzu1.a, gi.Gzu1.b, dx, dy, imgz1, ki);
Gzx1 = -1/j*shiftdim(dGzu1(1,:), 1); % Notice minus - this is because
Gzy1 = -1/j*shiftdim(dGzu1(2,:), 1); % differentiation by dx/dy in spatial
                                     % domain corresponds to multiplication
                                     % by -j*kx / -j*ky in spectral domain.

% Notice minus - because we need derivative with respect to zobs
[ Gzu2 dGzu2  ] = sumdcimg(gi.Gzu2.a, -gi.Gzu2.b, dx, dy, -imgz2, ki);
Gzx2 = -1/j*shiftdim(dGzu2(1,:), 1);
Gzy2 = -1/j*shiftdim(dGzu2(2,:), 1);

[ Gzu3 dGzu3  ] = sumdcimg(gi.Gzu3.a, gi.Gzu3.b, dx, dy, imgz3, ki);
Gzx3 = -1/j*shiftdim(dGzu3(1,:), 1);
Gzy3 = -1/j*shiftdim(dGzu3(2,:), 1);

% Notice minus - because we need derivative with respect to zobs
[ Gzu4 dGzu4  ] = sumdcimg(gi.Gzu4.a, -gi.Gzu4.b, dx, dy, -imgz4, ki);
Gzx4 = -1/j*shiftdim(dGzu4(1,:), 1);
Gzy4 = -1/j*shiftdim(dGzu4(2,:), 1);

Gvv = Gvv0+Gvv1+Gvv2+Gvv3+Gvv4;
Gzz = Gzz0+Gzz1+Gzz2+Gzz3+Gzz4;
Gzx =      Gzx1+Gzx2+Gzx3+Gzx4;
Gzy =      Gzy1+Gzy2+Gzy3+Gzy4;

