function fxy = integq_fxy(k, r, robs, qN, aux)
% function fxy = integq_fxy(k, r, robs, qN, aux)
% 
% Evaluates integral of div(fxy*exp(-j*k*R)/R) where fxy = f-z(f*z)
% (f without z component) using a triangle quadrature. Can not handle
% self-terms properly, used for the tests.
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
%  fxy =  d(fx*p)/dx + d(fy*p)/dy
%
% To evaluate the integral, we divide it into two parts. The first part
% includes the observation coordinate only, while the second depends
% on the free vertex position.
%

[ qdrx, qdry, qdrz ] = uncat(2, aux.qdr);
p = calcp(k, qdrx, qdry, qdrz);
dp = diffp(k, qdrx, qdry, qdrz);
dpdx = shiftdim(dp(1,:,:,:), 1);
dpdy = shiftdim(dp(2,:,:,:), 1);
dpdz = shiftdim(dp(3,:,:,:), 1);

[ qrx, qry, qrz ] = uncat(2, aux.qr);

%
% First part - depends on the observation coordinate only.
% f without Z component
divfxy = 2-sum(aux.n(:,:,1:2).^2, 3); % divergence of fxy
integrand = qrx.*dpdx+qry.*dpdy+p.*divfxy(:,:,ones(1,size(p, 3)));
qW = permute(aux.qW, [ 2 3 1 ]);
fxy_obs = permute(sum(integrand.*repmat(qW, N, 1), 3), [ 1 3 2 ]);

% Second part - includes the free vertex position.
intg_dpdx = sum(dpdx.*repmat(qW, N, 1), 3);
intg_dpdy = sum(dpdy.*repmat(qW, N, 1), 3);

[ rx, ry, rz ] = uncat(3, r);

fxy_r123 = rx.*intg_dpdx(:,ones(1,3))+ry.*intg_dpdy(:,ones(1,3));

% We are done! - just sum the two parts, and multiply by the areas
fxy = (repmat(fxy_obs, 1, 3) - fxy_r123).*aux.A(:,ones(1,3))*2;
