function ip = integq_p(k, r, robs, qN, aux)
% 
% Computes the same integral as computed by integ_p using a surface
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

qW = permute(aux.qW, [ 2 3 1 ]);
p = exp(-j*k*aux.qR)./aux.qR;
intgp = sum(p.*repmat(qW, N, 1), 3);
ip = intgp.*aux.A*2;
