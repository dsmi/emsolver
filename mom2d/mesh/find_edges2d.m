function e = find_edges2d(edges, verts, x, y, radius)
% e = find_edges2d(edges, verts, x, y, radius)
%
% Finds the mesh edges located within a particular radius around a given
% point.
% Inputs:
%    edges  - num_of_edges-by-2 matrix of the indices of the edge endpoint
%             vertices
%    verts  - num_of_verts-by-2 matrix of the coordinates of the vertices.
%    x, y   - center point of the lookup area
%    radius - radius of the lookup area
% Outputs:
%    e      - vector of the found edge indices

% Distance from the given point
r = sqrt(sum((verts - repmat([ x y ], size(verts, 1), 1)).^2,2));

% Distances for the edge points.
r1  = r(edges(:,1));
r2  = r(edges(:,2));

% Find the conforming edges.
e = find(abs(r1) <= radius & abs(r2) <= radius);
