function test_find_edges2d
% Test find_edges2d function, which finds the edges within a given radius from
% a given point.
%

% The testing geometry - a cicrle
[ edges, verts ] = mkcir2d(1, 10);

% Move it!
verts = verts + repmat([ 3 5 ], size(verts, 1), 1);

assertEquals(0, length(find_edges2d(edges, verts, 3, 5, 0.9)));
assertEquals(size(edges, 1), length(find_edges2d(edges, verts, 3, 5, 1.1)));
