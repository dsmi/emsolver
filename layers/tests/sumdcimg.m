function [ g dg d2g ] = sumdcimg(a, b, x, y, z, ki)
% [ g dg d2g ] = sumdcimg(a, b, x, y, z, ki)
%
% Used by test_sommerf - obtains the spatial domain value of a certain
% component of the layered Green's function by applying the Sommerfeld
% identity to the exponentials approximating the spectral domain value.
% In dcim, the spectral domain layered Greens function is approximated
% with a finite sum of exponentials, and then the approximating exponentials
% are transformed to the spatial domain using the Sommerfeld identity.
% Also returns the derivative of approximating exponentials with respect
% to x, y and z, it is used for quantities which have a linear dependence
% on kx and ky. For such quantities the following approach is used.
% kx/ky multiplier is dropped, and the function is appoximated with
% exponentials. Then, the exponentials (complex images) get converted
% to the spatial domain using the Sommerfeld identity; after that,
% the exponent is differentiated by x/y, and that finally gives the spatial
% domain value, because the differentiation by x/y in spatial domain
% corresponds to multiplication by j*kx/j*ky in spectral domain according to
% the Fourier transform properties.
% Inputs:
%   a       - multipliers of the approximating complex images.
%   x, y, z - Coordinates of the complex images relative to
%             the observation position.
%   ki      - wavenumber in the source/observation layer.
% Output:
%   d       - the resulting spatial domain value.
%   dg      - derivative of G with respect to x, y and z
%             size(dg) = [ 3 size(g) ]
%   d2g     - second derivative of G
%             size(d2g) = [ 3 3 size(g) ]
%

xr = repmat(shiftdim(x, -1), length(a), 1);
yr = repmat(shiftdim(y, -1), length(a), 1);
zr = repmat(shiftdim(z, -1), length(a), 1);
ar = repmat(a(:), [ 1 size(x) ]);
br = repmat(b(:), [ 1 size(x) ]);
zb = zr + br;
p = calcp(ki,xr,yr,zb);
g = 1/(2*pi)*shiftdim(sum(ar.*p, 1), -1);
dp = diffp(ki,xr,yr,zb);
dg = 1/(2*pi)*sum(repmat(shiftdim(ar,-1),3,1).*dp, 2);
d2p = diff2p(ki,xr,yr,zb);
d2g = 1/(2*pi)*sum(repmat(shiftdim(ar,-2),3,3).*d2p, 3);
