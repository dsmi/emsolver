function cfxy = integq_cfxy(k, r, robs, qN, aux)
% function cfxy = integq_cfxy(k, r, robs, qN, aux)
% 
% Evaluates integral of curl(z(fx*dg/dx+fy*dg/dy)) where fx and fy are
% the x and y componets of the basis function f = r - r(l) using a surface
% quadrature. Can not handle self-terms properly, used for the tests.
%

% Number of the quadrature points to use
if nargin < 4
	qN = 3;
end

% Auxiliaries, if not supplied by caller
if nargin < 5
	aux = integq_aux(r, robs, qN);
end

% Number of source-observation pairs.
N = size(r,1);

% We compute the integral of:
%  cfxy_x =  fx*d2g/dydx+fy*d2g/dydy
%  cfxy_y = -(fx*d2g/dxdx+fy*d2g/dxdy)
%
% To evaluate the integral, we divide it into two parts. The first part
% includes the observation coordinate only, while the second depends
% on the free vertex position.
%

[ qdrx, qdry, qdrz ] = uncat(2, aux.qdr);
d2p = diff2p(k, qdrx, qdry, qdrz);
d2p_dxdx = shiftdim(d2p(1,1,:,:,:), 2);
d2p_dxdy = shiftdim(d2p(1,2,:,:,:), 2);
d2p_dydx = d2p_dxdy;
d2p_dydy = shiftdim(d2p(2,2,:,:,:), 2);

[ qrx, qry, qrz ] = uncat(2, aux.qr);

%
% First part - depends on the observation coordinate only.
ix = qrx.*d2p_dydx+qry.*d2p_dydy;
iy = -(qrx.*d2p_dxdx+qry.*d2p_dxdy);
integrand = cat(2, ix, iy, 0*ix);
qW = permute(aux.qW, [ 2 3 1 ]);
cfxy_obs = permute(sum(integrand.*repmat(qW, N, 3), 3), [ 1 3 2 ]);

% Second part - includes the free vertex position.
intg_d2p_dxdx = sum(d2p_dxdx.*repmat(qW, N, 1), 3);
intg_d2p_dxdy = sum(d2p_dxdy.*repmat(qW, N, 1), 3);
intg_d2p_dydx = intg_d2p_dxdy;
intg_d2p_dydy = sum(d2p_dydy.*repmat(qW, N, 1), 3);

[ rx, ry, rz ] = uncat(3, r);

cfxy_r123_x = rx.*intg_d2p_dydx(:,ones(1,3))+ry.*intg_d2p_dydy(:,ones(1,3));
cfxy_r123_y = -(rx.*intg_d2p_dxdx(:,ones(1,3))+ry.*intg_d2p_dxdy(:,ones(1,3)));

cfxy_r123 = cat(3, cfxy_r123_x, cfxy_r123_y, 0*cfxy_r123_x);

% We are done! - just sum the two parts, and multiply by the areas
cfxy = (repmat(cfxy_obs, 1, 3) - cfxy_r123).*aux.A(:,ones(1,3),ones(1,3))*2;
