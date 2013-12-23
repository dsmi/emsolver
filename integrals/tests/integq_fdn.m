function fdn = integq_fdn(k, r, robs, qN, aux)
% function fdn = integq_fdn(k, r, robs, qN, aux)
% 
% Evaluates integral of f*d(exp(-j*k*R)/R)/dn using a triangle quadrature.
% Can not handle self-terms properly, used for the tests.
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
%  fdn =  f*dot(n, exp(-j*k*R)/R)
%
% To evaluate the integral, we divide it into two parts. The first part
% includes the observation coordinate only, while the second depends
% on the free vertex position.
%

[ qdrx, qdry, qdrz ] = uncat(2, aux.qdr);
dp = diffp(k, qdrx, qdry, qdrz);
dpdn = sum(permute(dp, [ 2 3 1 4 ]).*repmat(aux.n, [ 1 1 1 size(dp, 4) ]), 3);

%
% First part - depends on the observation coordinate only.
integrand = permute(aux.qr, [ 1 4 2 3 ]).*repmat(dpdn, [ 1 1 3 ]);
qW = shiftdim(aux.qW, -3);
fdn_obs = sum(integrand.*repmat(qW, [ N 1 3 ]), 4);

% Second part - includes the free vertex position.
intg_dpdn = sum(dpdn.*repmat(qW, N, 1), 4);

fdn_r123 = r.*intg_dpdn(:,ones(1,3),ones(1,3));

% We are done! - just sum the two parts, and multiply by the areas
fdn = (repmat(fdn_obs, 1, 3) - fdn_r123).*repmat(aux.A, [ 1 3 3 ])*2;
