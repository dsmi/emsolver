function test_shapes
%
% The mesh may consist of more than one separate shapes, and one of the tasks
% if the init_mesh function is to build lists of faces and edges which belong
% to a particular shape. Test this.
%

% first, try a mesh consisting of one shape
[ t, x, y, z ] = mkring(1, 0.1, 6, 6);

mesh = init_mesh(t, x, y, z);

assertEquals((1:size(t, 1))', mesh.shape_tris{1});
assertEquals((1:size(mesh.edges, 1))', mesh.shape_edges{1});

% next, the geometry composed of three shapes
[ t1, x1, y1, z1 ] = mkring(1, 0.1, 6, 6);
[ t2, x2, y2, z2 ] = mkring(0.1, 0.01, 8, 8);
[ t3, x3, y3, z3 ] = mkring(0.01, 0.001, 10, 10);
[ t, x, y, z ] = joinmeshes( { t1, t2, t3 }, { x1, x2, x3 }, { y1, y2, y3 }, { z1, z2, z3 } );

mesh = init_mesh(t, x, y, z);

st1 = (1:size(t1, 1))';
st2 = length(t1)+(1:size(t2, 1))';
st3 = length(t1)+length(t2)+(1:size(t3, 1))';
assertEquals(st1, mesh.shape_tris{1});
assertEquals(st2, mesh.shape_tris{2});
assertEquals(st3, mesh.shape_tris{3});

% Manually find edges which belong to each primitive
trie1 = mesh.tri_edges(st1,:);
edges1 = unique(trie1(:));
trie2 = mesh.tri_edges(st2,:);
edges2 = unique(trie2(:));
trie3 = mesh.tri_edges(st3,:);
edges3 = unique(trie3(:));

assertEquals(edges1, mesh.shape_edges{1});
assertEquals(edges2, mesh.shape_edges{2});
assertEquals(edges3, mesh.shape_edges{3});
