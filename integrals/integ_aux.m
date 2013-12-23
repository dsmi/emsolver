function iaux = integ_aux(r, robs, qN)
% function iaux = integ_aux(r, robs, qN)
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

% Normalized normals. N-by-1-by-3 array.
n = nn ./ nl(:,:,ones(1,3));

% triangle area, temporary, for the tests
nnn = cross(edges(:,1,:),edges(:,2,:),3);
A = sqrt(sum(nnn.^2, 3))./2;

% Signed distance from observation point to the plane of triangle.
% Column vector.
%h = dot(r(:,1,:), n, 3) - dot(robs, permute(n, [ 1 3 2 ]), 2);
h = sum(r(:,1,:).*n, 3) - sum(robs.*permute(n, [ 1 3 2 ]), 2);
abs_h = sqrt(h.*h); % abs() can not handle complex images

% normal derivative of abs_h
d_abs_h_dn = h*0;
nz = find(abs(abs_h) > 1e-15);
d_abs_h_dn(nz) = h(nz)./abs_h(nz);

% Projection of the observation point to the plane of triangle. N-by-1-by-3
% array, first index is the triangle index, third index is X-Y-Z.
rho = permute(robs, [ 1 3 2 ]) + n.*h(:,:,ones(1,3));

% n repeated so it has the same dimensions as r, r1, r2
n_r = repmat(n, 1, 3);

% Edge normal, cross product of the edge tangent and triangle normal.
% First index is the triangle index, second index is the local index of
% the egde in triangle, third index is X-Y-Z.
nu = cross(t, n_r, 3);

% Quadrature point positions and weights for the integrals over edges.
[qX,qW] = GLNodeWt(qN);
%[qX,qW] = GLTable(qN);

qW = repmat(permute(qW, [ 2 3 1 ]), [ N 3 ]);

% Quadrature point positions. N-by-3-by-3-by-qN array. First index is the
% triangle index, second one is the local index of the egde in triangle,
% third index is X-Y-Z, fourth is the quadrature point index.
qr = repmat(r1,[ 1 1 1 qN ]) + repmat(t,[ 1 1 1 qN ]) ...
		.* repmat(s,[ 1 1 3 qN ]) ...
		  .* repmat(permute((qX*0.5 + 0.5),[ 2 3 4 1 ]), [ N 3 3 1 ]);

% Distance from the observation point to quadrature point
qR = sqrt(sum((qr - repmat(permute(robs, [ 1 3 2 ]),[ 1 3 1 qN ])).^2,3));

% Permute it so the third index is a quadrature point index.
qR = permute(qR, [ 1 2 4 3 ]);

% abs_h matrix repeated so it has the same dimensions as qR
abs_h_qr = repmat(abs_h, [ 1 3 qN ]);

% h matrix repeated - N-by-3 array, N is the triangle index
h_r = h(:,ones(1,3));

% rho repeated so it has the same dimensions as r, r1, r2
rho_r = repmat(rho, 1, 3);

% Vector from 'rho' point to the quadrature point.
q_rho = qr - repmat(rho_r, [ 1 1 1 qN ]);

% Squared distance from 'rho' point to the quadrature point, N-by-3-by-qN
% array.
q_rho_2 = permute(sum(q_rho.^2,3), [1 2 4 3]);

r2_minus_rho = r2 - rho_r;
%rho2line = dot(r2_minus_rho, nu, 3);
rho2line = sum(r2_minus_rho.*nu, 3);

% Point-in-triangle test
%edge_x_rho = cross(edges,-r2_minus_rho,3);
%x_dot_n = sum(edge_x_rho.*repmat(nn, 1, 3), 3);
%cin = sum(x_dot_n >= 0, 2) > 2.9;

% The cin determines if the de-singularization term is needed, which is
% the case if rho is within or near the triangle.
tc = sum(r, 2)/3;              % triangle centers
rho_r2 = sum((rho-tc).^2, 3);  % squared distance from rho to the center
cin = rho_r2 < (max(s,[],2)*10).^2;

iaux.t = t;
iaux.n = n;
iaux.A = A; % temporary, for the tests
iaux.q_rho = q_rho; % temporary, for the tests
iaux.rho = rho; % temporary, for the tests
iaux.nu = nu;
iaux.s = s;
iaux.h = h;
iaux.rho2line = rho2line;
iaux.h_r = h_r;
iaux.n_r = n_r;
iaux.rho_r = rho_r;
iaux.abs_h = abs_h;
iaux.d_abs_h_dn = d_abs_h_dn;
iaux.qr = qr;
iaux.qR = qR;
iaux.qW = qW;
iaux.q_rho_2 = q_rho_2;
iaux.cin = cin;
