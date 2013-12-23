function test_find_eol2d
% Test find_edges2d function, which finds edges on the line
%

% The testing geometry - a cicrle
[ edges, verts ] = mkrect2d(1, 2, 5, 5);

e = find_eol2d(edges, verts, 1, 0, 0.5);
vi = edges(e,:); 
assertEquals( 0.5, verts(vi(:), 1), 1e-15);

e = find_eol2d(edges, verts, 1, 0, -0.5);
vi = edges(e,:); 
assertEquals( -0.5, verts(vi(:), 1), 1e-15);

e = find_eol2d(edges, verts, 0, 1, 1);
vi = edges(e,:); 
assertEquals( 1, verts(vi(:), 2), 1e-15);

e = find_eol2d(edges, verts, 0, 1, -1);
vi = edges(e,:); 
assertEquals( -1, verts(vi(:), 2), 1e-15);
