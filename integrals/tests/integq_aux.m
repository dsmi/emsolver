function iaux = integq_aux(r, robs, qN)
% function iaux = integq_aux(r, robs)
% 
% Calculates auxiliary geometry dependent values used by integration
% routines.
%
% Inputs:
%   r    - vertices of the triangles, array of size N-by-3-by-3, first index
%          is the index of the source-observation pair, second index is the
%          triangle vertex index, third is X-Y-Z.
%   robs - observation points, array of size N-by-3, first index is 
%          the index of the source-observation pair, second index is X-Y-Z.
%   qN   - The target surface integral is evaluated by reduction to the
%          line integrals over edges. qN is the number of quadrature points
%          to use over each edge. This parameter is optional, the default
%          value is 3.
%
% Outputs:
%   iaux - Structure with the following fields:
%      
%

% Number of quadrature points for the line integrals over edges.
if nargin < 3
	qN = 64;
end

% Number of source-observation pairs.
N = size(r,1);

% Edge endpoints. First edge is one opposite to the first vertex (bounded by
% vertices 2 and 3), second is one between vertices 3 and 1, third is
% between 1 and 2. N-by-3-by-3 array. First index is the local index of the
% egde in triangle, second one is the triangle index, third index is X-Y-Z.
r1 = zeros(N, 3, 3);
r1(:,1,:) = r(:,2,:);
r1(:,2,:) = r(:,3,:);
r1(:,3,:) = r(:,1,:);
r2 = zeros(N, 3, 3);
r2(:,1,:) = r(:,3,:);
r2(:,2,:) = r(:,1,:);
r2(:,3,:) = r(:,2,:);

% Edge vectors. N-by-3-by-3 array. First index is the triangle index, second
% index is the local index of the egde in triangle, third index is X-Y-Z.
edges = r2 - r1;

% Edge lengths. N-by-3 array. First index is the triangle index, second one
% is the local index of the egde in triangle.
s = sqrt(sum(edges.^2,3));

% Edge tangentials - normalized edge vectors.
t = edges ./ s(:,:,ones(1,3));

% Triangle normals, non-normalized yet. N-by-1-by-3 array, first index is
% the triangle index, third index is X-Y-Z.
nn = cross(t(:,1,:),t(:,2,:),3);

% Length of the normals, column vector.
nl = sqrt(sum(nn.^2, 3));

% Normalized normals. N-by-3 array.
n = nn ./ nl(:,:,ones(1,3));

% triangle area, temporary, for the tests
nnn = cross(edges(:,1,:),edges(:,2,:),3);
A = sqrt(sum(nnn.^2, 3))./2;

% Signed distance from observation point to the plane of triangle.
% Column vector.
%h = dot(r(:,1,:), n, 3) - dot(robs, permute(n, [ 1 3 2 ]), 2);
h = sum(r(:,1,:).*n, 3) - sum(robs.*permute(n, [ 1 3 2 ]), 2);
abs_h = sqrt(h.*h); % abs() can not handle complex images

% Projection of the observation point to the plane of triangle. N-by-1-by-3
% array, first index is the triangle index, third index is X-Y-Z.
rho = permute(robs, [ 1 3 2 ]) + n.*h(:,:,ones(1,3));

% n repeated so it has the same dimensions as r, r1, r2
n_r = repmat(n, 1, 3);

% Edge normal, cross product of the edge tangent and triangle normal.
% First index is the triangle index, second index is the local index of
% the egde in triangle, third index is X-Y-Z.
nu = cross(t, n_r, 3);


qNN = qN.*qN;
[ qA qW ] = simplexquad(qN, 2);
qA = [ qA 1-sum(qA,2) ]; % barycentric coordinates
qA = repmat(permute(qA, [ 3 2 4 1 ]), [ N 1 3 ]);

% Quadrature point positions. N-by-3-by-qN array. First index is the
% triangle index, second index is X-Y-Z, third one is the quadrature
% point index.
qr = permute(sum(repmat(r, [ 1 1 1 qNN ]).*qA, 2), [ 1 3 4 2 ]);

% Vector from observation point to quadrature point
qdr = qr - repmat(robs, [ 1 1 qNN ]);

% Distance from the observation point to quadrature point
qR = sqrt(sum(qdr.^2,2));

% Permute it so the second index is a quadrature point index.
%qR = permute(qR, [ 1 3 2 ]);	

iaux.n = n;
iaux.A = A; % temporary, for the tests
iaux.nu = nu;
iaux.s = s;
iaux.h = h;
iaux.qr = qr;
iaux.qdr = qdr;
iaux.qR = qR;
iaux.qW = qW;
