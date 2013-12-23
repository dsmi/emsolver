function test_loopstar2
% 
% One more test of the loopstar basis. This time it is tested with complex
% shapes which have holes, which means that the loopstar have to find
% additional basis cycles to augment the vertex-anchored cycles.
% To check if the resulting basis is correct we test: (1) if the transform
% matrix is full-rank, which means that the basis is complete (2) if the
% divergence of the loops is zero.
%

% The source primitive - three loops, different sizes and resolutions
[ t1, x1, y1, z1 ] = mkring(1, 0.5, 3, 3);
[ t2, x2, y2, z2 ] = mkring(0.1, 0.05, 5, 5);
[ t3, x3, y3, z3 ] = mkring(0.01, 0.005, 5, 5);

[ tri, x, y, z ] = joinmeshes({ t1 t2 t3 }, { x1 x2 x3 }, { y1 y2 y3 }, { z1 z2 z3 });

mesh = init_mesh(tri, x, y, z);

[ IL, IS ] = mkloopstar(mesh);
T = [ IL IS ];

% Check if the transform matrix is full rank, which is the case only if
% the basis is complete.
assertEquals(size(T, 1), rank(T))

% The arbitrary loop basis function coefficients.
X = (1:size(IL, 2))';

X = IL*X;
D = mkdivmat(mesh);
face_div = D*X;

% make sure the divergence of the loops is zero
assertEquals(0, face_div, 1e-9)
