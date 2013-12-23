function d2p = diff2p(k, x, y, z)
% d2p = diff2p(k, x, y, z)
%
% Calculates second derivative of fundamental solution of
% the Helmholtz equation.
% size(d2p) = [ 3 3 size(x) ]
%

r = calcr(x,y,z);
p = calcp(k,x,y,z);
dr = diffr(x,y,z);
d2r = diff2r(x,y,z);
dp = diffp(k,x,y,z);

% dr/dx_j
drj = repmat(shiftdim(dr, -1), 3, 1);

% dr/dx_k, dp/dx_k, 
pvec = 1:(length(size(x))+2);
pvec(1:2) = [ 2 1 ];
drk = repmat(permute(shiftdim(dr, -1), pvec), 1, 3);
dpk = repmat(permute(shiftdim(dp, -1), pvec), 1, 3);

invr = repmat(shiftdim(1./r, -2), 3, 3);
invr2 = repmat(shiftdim(1./(r.*r), -2), 3, 3);
ppp = repmat(shiftdim(p, -2), 3, 3);

d2p = invr2.*ppp.*drk.*drj-(j*k+invr).*(dpk.*drj+ppp.*d2r);
