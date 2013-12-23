function cfxy = integ_cfxy(k, r, robs, qN, aux)
% function cfxy = integ_cfxy(k, r, robs, qN, aux)
% 
% Calculates integrals of the 
%   curl(z(fx*dg/dx+fy*dg/dy)), where
%    fx and fy are the x and y components of f = r-r(l)
%    g = exp(-j*k*R)/R
%    curl operates only on g
%
% Inputs:
%   k    - wavenumber, scalar
%   r    - vertices of the triangles, array of size N-by-3-by-3, first index
%          is the index of the source-observation pair, second index is the
%          triangle vertex index, third is X-Y-Z.
%   robs - observation points, array of size N-by-3, first index is 
%          the index of the source-observation pair, second index is X-Y-Z.
%   qN   - The target surface integral is evaluated by reduction to the
%          line integrals over edges. qN is the number of quadrature points
%          to use over each edge. This parameter is optional, the default
%          value is 3.
%   aux  - structure carrying auxiliary values, as returned by integ_aux.
%          When this function is called multiple times with the same r
%          and robs but with different k, it may be reasonable from
%          the performance point of view to compute the auxiliary values
%          once and use them for all the calls.
%          This parameter is optional, if aux is not supplied by the caller
%          integ_aux is called.
%   
% Outputs:
%   cfxy - the resulting integrals. Three-dimensional array of size
%          N-by-3-by-3, first index is the source-observation pair index,
%          second index is the local index of the free vertex,
%          third is X-Y-Z. Z is always zero!
%

% Number of the quadrature points to use
if nargin < 4
	qN = 3;
end

% Auxiliaries, if not supplied by caller
if nargin < 5
	aux = integ_aux(r, robs, qN);
end

% We compute the integral of:
%  cfxy_x =  fx*d2g/dydx+fy*d2g/dydy
%  cfxy_y = -(fx*d2g/dxdx+fy*d2g/dxdy)
%
% First, we compute the integrals of fx*grad(dp/dx), fy*grad(dp/dy)
%  where fx = x*dot(f,x), x is the basis vector, and p = exp(-j*k*R)/R.
% It can be computed as:
%  fx*grad(dp/dx) = grad(fx*dp/dx) - dp/dx*grad(fx)
% Next:
%  grad(fx*dp/dx) = grad_s(fx*dp/dx) + d/dn(fx*dp/dx), grad_s is the surface
%                    gradient.
%

% p = exp(-j*k*R)./R
p = exp(-j*k*aux.qR)./aux.qR;

% derivative with respect to x if multiplied by x and with respect to y if
% multiplied by y
dpd = -p.*(1+j*k*aux.qR)./(aux.qR.*aux.qR);

[ qrx qry qrz ] = uncat(4, permute(aux.qr, [ 1 2 4 3 ]));

qdr = aux.qr - repmat(permute(robs, [ 1 3 2 ]),[ 1 3 1 qN ]);
[ qdrx qdry qdrz ] = uncat(4, permute(qdr, [ 1 2 4 3 ]));

s2 = aux.s./2;

% Surface gradient of d(exp(-j*k*R)/R)/dx and d(exp(-j*k*R)/R)/dy
gs_dpdx = sum(aux.nu.*repmat(sum(dpd.*qdrx.*aux.qW,3).*s2, [ 1 1 3 ]), 2);
gs_dpdy = sum(aux.nu.*repmat(sum(dpd.*qdry.*aux.qW,3).*s2, [ 1 1 3 ]), 2);

% This is the part of the integral of grad_s(f(x/y)*dp/d(x/y)) which
% depends on the free vertex
gs_fxdpdx_r123 = repmat(r(:,:,1), [ 1 1 3 ]).*repmat(gs_dpdx, [ 1 3 ]);
gs_fydpdx_r123 = repmat(r(:,:,2), [ 1 1 3 ]).*repmat(gs_dpdx, [ 1 3 ]);
gs_fxdpdy_r123 = repmat(r(:,:,1), [ 1 1 3 ]).*repmat(gs_dpdy, [ 1 3 ]);
gs_fydpdy_r123 = repmat(r(:,:,2), [ 1 1 3 ]).*repmat(gs_dpdy, [ 1 3 ]);

% This is the part of the integral of grad_s(f(x/y)*dp/d(x/y))
% which depends on the observation point only
gs_fxdpdx_r = sum(aux.nu.*repmat(sum(dpd.*qrx.*qdrx.*aux.qW,3).*s2, [ 1 1 3 ]), 2);
gs_fydpdx_r = sum(aux.nu.*repmat(sum(dpd.*qry.*qdrx.*aux.qW,3).*s2, [ 1 1 3 ]), 2);
gs_fxdpdy_r = sum(aux.nu.*repmat(sum(dpd.*qrx.*qdry.*aux.qW,3).*s2, [ 1 1 3 ]), 2);
gs_fydpdy_r = sum(aux.nu.*repmat(sum(dpd.*qry.*qdry.*aux.qW,3).*s2, [ 1 1 3 ]), 2);

% This is the integral of grad_s(f(x/y)*dp/d(x/y))
gs_fxdpdx = gs_fxdpdx_r(:,ones(1,3),:) - gs_fxdpdx_r123;
gs_fydpdx = gs_fydpdx_r(:,ones(1,3),:) - gs_fydpdx_r123;
gs_fxdpdy = gs_fxdpdy_r(:,ones(1,3),:) - gs_fxdpdy_r123;
gs_fydpdy = gs_fydpdy_r(:,ones(1,3),:) - gs_fydpdy_r123;

% Integrals of d/dn(f(x/y)*dp/d(x/y))
fdndpdx = integ_fdndv(k, r, robs, [ 1 0 0 ], qN, aux);
fdndpdy = integ_fdndv(k, r, robs, [ 0 1 0 ], qN, aux);
ddn_fxdpdx = fdndpdx(:,:,[ 1 1 1 ]).*aux.n_r;
ddn_fydpdx = fdndpdx(:,:,[ 2 2 2 ]).*aux.n_r;
ddn_fxdpdy = fdndpdy(:,:,[ 1 1 1 ]).*aux.n_r;
ddn_fydpdy = fdndpdy(:,:,[ 2 2 2 ]).*aux.n_r;

% Gradients of fx, fy and fx
dfx = repmat(cat(3, 1, 0, 0), size(r, 1), 1) - aux.n(:,:,[ 1 1 1 ]).*aux.n;
dfy = repmat(cat(3, 0, 1, 0), size(r, 1), 1) - aux.n(:,:,[ 2 2 2 ]).*aux.n;

% Integrals of dp/d(x/y)*grad(f(x/y))
dp = integ_dp(k, r, robs, qN, aux);
dfx_dpdx = dp(:,:,[ 1 1 1 ]).*dfx;
dfy_dpdx = dp(:,:,[ 1 1 1 ]).*dfy;
dfx_dpdy = dp(:,:,[ 2 2 2 ]).*dfx;
dfy_dpdy = dp(:,:,[ 2 2 2 ]).*dfy;

fx_gdpdx = gs_fxdpdx - repmat(dfx_dpdx, 1, 3) + ddn_fxdpdx;
fy_gdpdx = gs_fydpdx - repmat(dfy_dpdx, 1, 3) + ddn_fydpdx;
fx_gdpdy = gs_fxdpdy - repmat(dfx_dpdy, 1, 3) + ddn_fxdpdy;
fy_gdpdy = gs_fydpdy - repmat(dfy_dpdy, 1, 3) + ddn_fydpdy;

% Finally, given the integrals of fx*grad(dp/dx), fy*grad(dp/dx),
% fx*grad(dp/dy) and fx*grad(dp/dy) compute curl(z(fx*dg/dx+fy*dg/dy))
cfxy = 0*fx_gdpdy;
cfxy(:,:,1) = fx_gdpdx(:,:,2) + fy_gdpdy(:,:,2); % fx*d2g/dydx+fy*d2g/dydy
cfxy(:,:,2) = -(fx_gdpdx(:,:,1) + fy_gdpdy(:,:,1)); % -(fx*d2g/dxdx+fy*d2g/dxdy)
