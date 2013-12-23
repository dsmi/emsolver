function test_ls_manyobj
% test_ls_manyobj : Test if mkloopstar works correctly if the geometry
% consists of a number of isolated objects.
%

% The test idea is to create the loop and star transform matrices for two
% primitives separately and compose the matrices, then create the matrices
% for the united geometry and compare the results.

% Parameters of the matter - free space
permittivity = eps0;
permeability = mu0;

% Angular frequency
freq = 1e5;

% Wavenumber used when evaluating integrals
k = freq * sqrt(permeability * permittivity);

% The first primitive is 1x1x1 box.
[ tri1, x1, y1, z1 ] = mkbox(1, 1, 1, 1, 1, 1);

mesh1 = init_mesh(tri1, x1, y1, z1);

[ IL1, IS1 ] = mkloopstar(mesh1);

% The second primitive - cylinder.
[ tri2, x2, y2, z2 ] = mkpole(1, 0.5, 5, 5, 1);
[ x2, y2, z2 ] = move(x2, y2, z2, 2, 0, 0);

mesh2 = init_mesh(tri2, x2, y2, z2);

[ IL2, IS2 ] = mkloopstar(mesh2);

% Now create a mesh which is union of the two primitives
tri = [ tri1; tri2 + length(x1) ];
x = [ x1 x2 ];
y = [ y1 y2 ];
z = [ z1 z2 ];

mesh = init_mesh(tri, x, y, z);

[ IL, IS ] = mkloopstar(mesh);

% The test matrices are created by joining matrices created for the test
% shapes separately.
test_IL = [           IL1           zeros(size(IL1,1), size(IL2,2)) ; ...
         zeros(size(IL2,1), size(IL1,2))          IL2    ];
test_IS = [           IS1           zeros(size(IS1,1), size(IS2,2)) ; ...
         zeros(size(IS2,1), size(IS1,2))          IS2    ];

assertEquals(test_IL, IL);
assertEquals(test_IS, IS);
