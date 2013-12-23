function fxg = integq_tfxg(k, r, robs, M, qN)
% function fxg = integq_tfxg(k, r, robs, M, qN)
% 
% Calculates integrals of the 
%   M*f X grad(exp(-j*k*R)/R),
% where f is the vector between free vertex and integration point: f = r-r(l),
% X means cross product, R is the distance between observation and integration
% points, and M is the caller-provided transform matrix, which is applied to f.
% The integral is evaluated using a surface quadrature. The function can not
% handle self-terms properly and it is used for the tests.
%

% Number of the quadrature points to use
if nargin < 5
	qN = 3;
end

% Auxiliaries, if not supplied by caller
if nargin < 6
	aux = integq_aux(r, robs, qN);
end

% Number of source-observation pairs.
N = size(r,1);

% We compute the integral of:
%  fxg_x = M(2,2)*fy*dgdz-M(3,3)*fz*dgdy
%  fxg_y = -M(1,1)*fx*dgdz+M(3,3)*fz*dgdx
%  fxg_z = M(1,1)*fx*dgdy-M(2,2)*fy*dgdx
%
% To evaluate the integral, we divide it into two parts. The first part
% includes the observation coordinate only, while the second depends
% on the free vertex position.
%

[ qdrx, qdry, qdrz ] = uncat(2, aux.qdr);
dp = diffp(k, qdrx, qdry, qdrz);
dpdx = shiftdim(dp(1,:,:,:), 1);
dpdy = shiftdim(dp(2,:,:,:), 1);
dpdz = shiftdim(dp(3,:,:,:), 1);

[ qrx, qry, qrz ] = uncat(2, aux.qr);

%
% First part - depends on the observation coordinate only.
ix = M(2,2)*qry.*dpdz-M(3,3)*qrz.*dpdy;
iy = -M(1,1)*qrx.*dpdz+M(3,3)*qrz.*dpdx;
iz = M(1,1)*qrx.*dpdy-M(2,2)*qry.*dpdx;
integrand = cat(2, ix, iy, iz);
qW = permute(aux.qW, [ 2 3 1 ]);
fxg_obs = permute(sum(integrand.*repmat(qW, N, 3), 3), [ 1 3 2 ]);

% Second part - includes the free vertex position.
intg_dpdx = sum(dpdx.*repmat(qW, N, 1), 3);
intg_dpdy = sum(dpdy.*repmat(qW, N, 1), 3);
intg_dpdz = sum(dpdz.*repmat(qW, N, 1), 3);

[ rx, ry, rz ] = uncat(3, r);

fxg_r123_x = M(2,2)*ry.*intg_dpdz(:,ones(1,3))-M(3,3)*rz.*intg_dpdy(:,ones(1,3));
fxg_r123_y = -M(1,1)*rx.*intg_dpdz(:,ones(1,3))+M(3,3)*rz.*intg_dpdx(:,ones(1,3));
fxg_r123_z = M(1,1)*rx.*intg_dpdy(:,ones(1,3))-M(2,2)*ry.*intg_dpdx(:,ones(1,3));

fxg_r123 = cat(3, fxg_r123_x, fxg_r123_y, fxg_r123_z);

% We are done! - just sum the two parts, and multiply by the areas
fxg = (repmat(fxg_obs, 1, 3) - fxg_r123).*aux.A(:,ones(1,3),ones(1,3))*2;
