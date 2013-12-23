function e = find_eol2d(edges, verts, nx, ny, d)
% e = find_eol2d(edges, verts, nx, ny, d)
%
% Finds the mesh edges lying on the given line.
%
% Inputs:
%    edges     - num_of_edges-by-2 matrix of the indices of the edge endpoint
%                vertices
%    verts     - num_of_verts-by-2 matrix of the coordinates of the vertices.
%    nx, ny, d - coefficients of the line equation of the form:
%                   nx*x+ny*y+d=0
%                (nx, ny) is the normal vector, d is the distance from zero
%                measured in the direction pointed by the normal
% Outputs:
%    e         - vector of the found edge indices

% Allow non-normalized normals
nl = sqrt(nx*nx+ny*ny);
nx = nx/nl;
ny = ny/nl;

% Distance from the given point
r = abs(sum((verts.*repmat([ nx ny ], size(verts, 1), 1)),2)-d);

% Distances for the edge points.
r1  = r(edges(:,1));
r2  = r(edges(:,2));

% Find the conforming edges.
tol = 1e-15;
e = find(abs(r1) <= tol & abs(r2) <= tol);
