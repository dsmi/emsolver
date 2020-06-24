function triv = calc_triv( mesh, xe )
% triv = calc_triv( mesh, xe )
%
% Given the weights of the RWG basis functions calculates values of the
% vector field at the triangle centers

% Value of the basis functions at the triangle centers scaled
edge_rc_s = mesh.edge_rc .* repmat( xe, [ 1 2 3 ] );

% Edges of a triangle, signs, and 123 indices for xyz
te  = repmat(mesh.tri_edges, [ 1 1 3 ]);
tes = repmat(mesh.tri_edges_s, [ 1 1 3 ]);
v   = repmat(shiftdim(1:3,-1), size( mesh.tri, 1 ), 3 );

% n-by-1-by-3
triv3 = sum( edge_rc_s( sub2ind( size(edge_rc_s), te, tes, v ) ), 2 );

% n-by-3
triv = permute( triv3, [ 1 3 2 ] );
