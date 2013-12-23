function fp = integq_fp(k, r, robs, qN, aux)
% function fp = integq_fp(k, r, robs, qN, aux)
% 
% Computes the same intagral as computed by integ_fp using a surface
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
% f*p = f*exp(-j*k*R)/R), where f is the RWG basis function given by
% f = r - r(l), and r is the evaluation point coordinates and r(l)
% is the free vertex.
%
% To evaluate the integral, we divide it into two parts:
%   x1 = r*exp(-j*k*R)/R
%   x2 = r(l)*exp(-j*k*R)/R
% The first part includes the observation coordinate only, while the second
% depends on the free vertex position as well.
%
% First part - depends on the observation coordinate only.
qW = permute(aux.qW, [ 2 3 1 ]);
p = exp(-j*k*aux.qR)./aux.qR;
fp_obs = permute(sum(aux.qr.*repmat(p, 1, 3).*repmat(qW, N, 3), 3), [ 1 3 2 ]);

% Second part - includes the free vertex position.
intgp = sum(p.*repmat(qW, N, 1), 3);
fp_r123 = r.*repmat(intgp, [ 1 3 3 ]);

% We are done! - just sum the two parts, and multiply by the areas
fp = (repmat(fp_obs, 1, 3) - fp_r123).*aux.A(:,ones(1,3),ones(1,3))*2;
