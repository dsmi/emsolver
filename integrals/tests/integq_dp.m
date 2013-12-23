function dp = integq_dp(k, r, robs, qN, aux)
% function dp = integ_dp(k, r, robs, qN, aux)
% 
% Evaluates integral of grad(exp(-j*k*R)/R) using a surface quadrature.
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

[ qdrx, qdry, qdrz ] = uncat(2, aux.qdr);
qdp = diffp(k, qdrx, qdry, qdrz);
qdp = permute(qdp, [ 2 3 1 4 ]);

dp = sum(qdp.*repmat(shiftdim(aux.qW, -3), [ N 1 3 ]), 4);

% multiply by the areas
dp = dp.*aux.A(:,:,ones(1,3))*2;
