function v = intg_helmdl2d(k,rsrc,robs)
% v = intg_helmdl2d(k,rsrc,robs)
%
% Evaluates integral of the double layer potential over the boundary
% segments, where the Green's function is the one associated with
% the 2D Helmholtz equation, which is given by
%   -j/4*besselh(0,2,k*r) (= 1/(2*pi)*besselk(0,j*k*r)).
% The double layer potential is the normal derivative of the Green's
% function given by the following formula
%   g*(R) = -j*k/(2*pi)*K1(j*k*r)*dr/dR*n
% where r is distance, R is position vector, n is normal, dr/dR = R/r.
% Notice that this function identifies singular source-observation pairs
% (ones where source coincides with observation) by comparison of the
% segment endpoints which only allows to identify the cases where the
% source segment exactly matches the observation one, any other singular
% configurations are ignored and may cause undefined behavior.
%
% Params:
%  k      - wavenumber.
%  rsrc   - source segments, N-by-2-by-2 array, idx-XYZ-begin/end.
%  robs   - observation segments, N-by-2-by-2 array, idx-XYZ-begin/end.
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

% Outer normals. The outer polygon is CCW, the holes are CW, thus
% to get the outer normals we need to rotate the edge tangentials
% 90 degrees clockwise.
n = [ t(:,2) -t(:,1) ];

% Edge centers, N-by-2
obsc = sum(robs, 3)*0.5;

% Find singular pairs
sidx = find(all(all(abs(rsrc-robs)<1e-12, 3), 2));

% Quadrature point positions. N-by-2-by-qN array. First index is the
% edge index, second one is X-Y, third is the quadrature point index.
qr = repmat(rsrc(:,:,1), [ 1 1 qN ]) + repmat(edges, [ 1 1 qN ]) ...
            .* repmat(permute(qX*0.5+0.5, [ 2 3 1 ]), N, 2);

% Now it is a vector from the observation point to the quadrature point.
qr = qr - repmat(obsc,[ 1 1 qN ]);

% Distance from the observation point to quadrature point. N-by-1-by-qN.
qR = sqrt(sum(qr.^2,2));

dr = qr./repmat(qR, 1, 2); % Derivative of the distance.
dr_dot_n = permute(dot(dr, repmat(n, [ 1 1 qN ]), 2), [ 1 3 2 ]);

% Permute qR it so the second index is the quadrature point index.
qR = permute(qR, [ 1 3 2 ]);

% Calculate the quadratures.
integrand1 = -j*k/(2*pi)*besselk(1,j*k*qR).*dr_dot_n;
qW_2 = repmat(permute(qW, [ 2 1 ]), N, 1);
integral1 = sum(integrand1.*qW_2,2).*(s/2);

v = integral1;

% For singular elements, this integral is zero because it involves the
% normal derivative of the distance.
v(sidx) = 0;
