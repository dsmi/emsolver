function test_intg_laplace2d
%
% Test test_intg_lapsl2d and test_intg_lapsl2d routines, which compute
% single and double layer potential of the Green's function of
% Laplace equation. The functions are tested by comparison of
% results with a fine quadrature.
%

% The testing geometry - a cicrle
[ edges, verts ] = mkcir2d(1, 100);
rsrc = cat(3, verts(edges(:,1),:), verts(edges(:,2),:));
robs = repmat([ linspace(-1, 1, 100)' linspace(1, -1, 100)' ], [ 1 1 2 ]);

% Plus a singular edge test
rsrc = [ rsrc ; cat(3, [ 0 0 ], [ 1 1 ]) ];
robs = [ robs ; cat(3, [ 0 0 ], [ 1 1 ]) ];

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

% Number of quadrature points used.
qN=100;

[qX,qW] = GLNodeWt(qN);

% Number of edges
N = size(rsrc,1);

% Edge centers, N-by-2
obsc = sum(robs, 3)*0.5;

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

% Permute qR so the second index is the quadrature point index.
qR = permute(qR, [ 1 3 2 ]);

% Calculate quadrature for the single layer kernel.
integrand = -1/(2*pi)*log(qR);
qW_2 = repmat(permute(qW, [ 2 1 ]), N, 1);
sl_test = sum(integrand.*qW_2,2).*(s/2);
sl = intg_lapsl2d(rsrc, robs);
assertEquals(sl_test(1:100), sl(1:100), 1e-15);
assertEquals(sl_test(101), sl(101), 5e-3);

% Calculate quadrature for the double layer kernel.
integrand = -1./(2*pi*qR).*dr_dot_n;
dl_test = sum(integrand.*qW_2,2).*(s/2);
dl = intg_lapdl2d(rsrc, robs);
assertEquals(dl_test(1:100), dl(1:100), 1e-15);
assertEquals(dl_test(101), dl(101), 1e-8);
