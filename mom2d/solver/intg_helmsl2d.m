function v = intg_helmsl2d(k,rsrc,robs)
% v = intg_helmsl2d(k,rsrc,robs)
%
% Evaluates integral of the single layer potential over the boundary
% segments, where the Green's function is the one associated with
% the 2D Helmholtz equation, which is given by
% -j/4*besselh(0,2,k*r) (= 1/(2*pi)*besselk(0,j*k*r)).
%
% Params:
%  k      - wavenumber.
%  r1, r2 - edge endpoints, N-by-2 arrays.
%  robs   - observation points, N-by-2 array.
%  sidx   - indices of the singular pairs, i.e. pairs where observation
%           points lies on the edge. For such pairs, the Greens function
%           is replaced by small argument formula and the integral
%           is evaluated analytically.
%           

% Number of edges
N = size(rsrc,1);

% Number of quadrature points used.
qN=7;

% Quadrature point positions and weights for the integrals over edges.
%[qX,qW] = GLNodeWt(qN);
[qX,qW] = GLTable(qN);

% Edge vectors. N-by-2 array. First index is the edge index, second
% one is X-Y.
edges = rsrc(:,:,2)-rsrc(:,:,1);

% Edge lengths. Column vector of length N.
s = sqrt(sum(edges.^2,2));

% Edge tangentials - normalized edge vectors.
t = edges ./ s(:,ones(1,2));

% Edge centers, N-by-2
obsc = sum(robs, 3)*0.5;

% Find singular pairs
sidx = find(all(all(abs(rsrc-robs)<1e-12, 3), 2));

% Quadrature point positions. N-by-2-by-qN array. First index is the
% edge index, second one is X-Y, third is the quadrature point index.
qr = repmat(rsrc(:,:,1), [ 1 1 qN ]) + repmat(edges, [ 1 1 qN ]) ...
            .* repmat(permute(qX*0.5+0.5, [ 2 3 1 ]), N, 2);

% Distance from the observation point to quadrature point. N-by-1-by-qN.
qR = sqrt(sum((qr - repmat(obsc,[ 1 1 qN ])).^2,2));

% Permute it so the second index is the quadrature point index.
qR = permute(qR, [ 1 3 2 ]);

% Calculate the quadratures.
integrand1 = -j/4*besselh(0,2,k*qR);
qW_2 = repmat(permute(qW, [ 2 1 ]), N, 1);
integral1 = sum(integrand1.*qW_2,2).*(s/2);

v = integral1;

% Euler's constant
gamma = 1.781072418;

% For singular elements, use the small argument formula for Hankel
% function (H02(z) = 1-j*2/pi*log(gamma*z/2)) and inegrate it analytically
v(sidx) = -j/4*s(sidx).*(1-j*2/pi*log((gamma*k*s(sidx))/(4*exp(1))));
