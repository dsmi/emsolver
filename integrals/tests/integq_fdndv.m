function fdndv = integq_fdndv(k, r, robs, v, qN, aux)
% function fdndv = integq_fdndv(k, r, robs, v, qN, aux)
% 
% Evaluates integral of f*d/dn(d/dv(p)) - the same as computed by integ_fdndv
% using a surface quadrature.
% Can not handle self-terms properly, used for the tests.
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

%
% To evaluate the integral, we divide it into two parts. The first part
% includes the observation coordinate only, while the second depends
% on the free vertex position.
%

v = reshape(v, 3, 1);

[ qdrx, qdry, qdrz ] = uncat(2, aux.qdr);
d2p = diff2p(k, qdrx, qdry, qdrz);
d2p_dv = shiftdim(sum(d2p.*repmat(v, [ 1 3 N 1 qN*qN ]), 1), 1);
d2p_dndv = shiftdim(sum(d2p_dv.*repmat(permute(aux.n, [ 3 1 2 ]), [ 1 1 1 qN*qN ]), 1), 1);

%
% First part - depends on the observation coordinate only.
integrand = aux.qr.*repmat(d2p_dndv, 1, 3);
qW = permute(aux.qW, [ 2 3 1 ]);
fdndv_obs = permute(sum(integrand.*repmat(qW, N, 3), 3), [ 1 3 2 ]);

% Second part - includes the free vertex position.
fdndv_r123 = r.*repmat(sum(d2p_dndv.*repmat(qW, N, 1), 3), [ 1 3 3 ]);

% We are done! - just sum the two parts, and multiply by the areas
fdndv = (repmat(fdndv_obs, 1, 3) - fdndv_r123).*aux.A(:,ones(1,3),ones(1,3))*2;
