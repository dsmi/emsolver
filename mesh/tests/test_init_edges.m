function test_init_edges
% test_init_edges : Test if init_edges does the job correctly.
%

% The source primitive
[ tri, x, y, z ] = mkring(5, 1, 20, 20);
%[ tri, x, y, z ] = mkbox(1, 1, 1, 1, 1, 1);

[ edges, etris, freev, freevl, trie, tries ] = init_edges(tri);

% Check the number of edges first.
ntri = size(tri,1);
nedges = size(edges,1);
assertEquals(ntri*3/2, nedges);

% Triangle centers, row vectors.
cx = sum(x(tri),2)'/3;
cy = sum(y(tri),2)'/3;
cz = sum(z(tri),2)'/3;

% Edge vectors, used to compute normals
e12 = [ x(tri(:,2)) - x(tri(:,1)); ...
        y(tri(:,2)) - y(tri(:,1)); ...
        z(tri(:,2)) - z(tri(:,1)) ];
e23 = [ x(tri(:,3)) - x(tri(:,2)); ...
        y(tri(:,3)) - y(tri(:,2)); ...
        z(tri(:,3)) - z(tri(:,2)) ];

% Triangle normals, 3-by-N array
tn = cross(e12,e23,1);

% A basic test - there should be no repeated vertices in a triangle.
assertTrue(isempty(find(edges(:,1) == edges(:,2))));
assertTrue(isempty(find(freev(:,1) == edges(:,1))));
assertTrue(isempty(find(freev(:,2) == edges(:,1))));
assertTrue(isempty(find(freev(:,1) == edges(:,2))));
assertTrue(isempty(find(freev(:,2) == edges(:,2))));

% For each edge, compute centers and normals of positive and negative
% triangles using edges and freev.
pos_cx = sum([ x(edges) x(freev(:,1))' ],2)'/3;
pos_cy = sum([ y(edges) y(freev(:,1))' ],2)'/3;
pos_cz = sum([ z(edges) z(freev(:,1))' ],2)'/3;
neg_cx = sum([ x(edges) x(freev(:,2))' ],2)'/3;
neg_cy = sum([ y(edges) y(freev(:,2))' ],2)'/3;
neg_cz = sum([ z(edges) z(freev(:,2))' ],2)'/3;

pos_e12 = [ x(edges(:,2)) - x(edges(:,1)); ...
            y(edges(:,2)) - y(edges(:,1)); ...
            z(edges(:,2)) - z(edges(:,1)) ];
pos_e23 = [ x(freev(:,1)) - x(edges(:,2)); ...
            y(freev(:,1)) - y(edges(:,2)); ...
            z(freev(:,1)) - z(edges(:,2)) ];

pos_tn = cross(pos_e12,pos_e23,1);

neg_e12 = [ x(edges(:,1)) - x(edges(:,2)); ...
            y(edges(:,1)) - y(edges(:,2)); ...
            z(edges(:,1)) - z(edges(:,2)) ];
neg_e23 = [ x(freev(:,2)) - x(edges(:,1)); ...
            y(freev(:,2)) - y(edges(:,1)); ...
            z(freev(:,2)) - z(edges(:,1)) ];

neg_tn = cross(neg_e12,neg_e23,1);

% Validate the calculated values against the triangles referred
% by etris
tol = 1e-14;
assertEquals(pos_cx, cx(etris(:,1)), tol);
assertEquals(pos_cy, cy(etris(:,1)), tol);
assertEquals(pos_cz, cz(etris(:,1)), tol);
assertEquals(neg_cx, cx(etris(:,2)), tol);
assertEquals(neg_cy, cy(etris(:,2)), tol);
assertEquals(neg_cz, cz(etris(:,2)), tol);
assertEquals(pos_tn, tn(:,etris(:,1)), tol);
assertEquals(neg_tn, tn(:,etris(:,2)), tol);

% Test if freevl is valid based on freev
idx = sub2ind(size(tri),etris,freevl);
assertEquals(tri(idx),freev);

% Validate trie and tries
for iedge=1:3,
  edge_v = repmat([ 1 2 ], ntri, 1);
  neg_edges = find(tries(:,iedge) == 2);
  edge_v(neg_edges,1:2) = repmat([ 2 1 ], length(neg_edges), 1);
  idx = sub2ind(size(edges),repmat(trie(:,iedge),1,2),edge_v);
  edge_v1 = [ 2 3 1 ];
  edge_v2 = [ 3 1 2 ];
  tri_edges = [ tri(:,edge_v1(iedge)) tri(:,edge_v2(iedge)) ];
  assertEquals(edges(idx),tri_edges);
end
