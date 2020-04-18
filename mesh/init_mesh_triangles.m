function mesh = init_mesh_triangles(tri, x, y, z)
% mesh = init_mesh_triangles(tri, x, y, z)
%
% Initializes the mesh structure, which contains a number of auxiliary
% values needed by for example the integration routines.
% This one is used when the solver is use works with triangles only and
% does not need edges. See init_mesh for the mesh structure which defines
% the edges as well.
%
%  Inputs: 
%    tri   -   triangles, num_of_tri-by-3 array, three vertex indices for each
%              triangle. The triangles need to be  oriented CCW when looking
%              from outside the geometry.
%    x, y, z - coordinates of the mesh vertices.
%
%  Output is a structure with the following fields:
%    tri     - triangles, as specified by the caller
%    x, y, z - vertices, as specified by the caller
%    cx, cy, cz - triangle centers, column vector of length num_of_triangles.
%    rc        - vectors from a vertex of the triangle to the triangle center.
%              3-by-num_of_triangles-by-3 array; first index is the local
%              index of the vertex in triangle, second is the triangle index,
%              third index is X-Y-Z.
%    nx, ny, nz - triangle normals, normalized, column vector of length
%              num_of_triangles. Oriented according to the right-hand rule.
%    tri_a   - Triangle areas, column vector of length num_of_triangles.
%

mesh.tri = tri;
mesh.x = x;
mesh.y = y;
mesh.z = z;


% Number of the triangles
ntri = size(tri,1);

% Triangle vertices
tvx = x(tri);
tvy = y(tri);
tvz = z(tri);

% Triangle centers
mesh.cx = sum(tvx,2)/3;
mesh.cy = sum(tvy,2)/3;
mesh.cz = sum(tvz,2)/3;

%  Vectors from a vertex to the center
mesh.rc = cat(3, (repmat(mesh.cx, 1, 3) - tvx).', (repmat(mesh.cy, 1, 3) - tvy).', ...
            (repmat(mesh.cz, 1, 3) - tvz).');

% Compute the normals as the cross product of the edges
edge12 = [ tvx(:,2) - tvx(:,1), tvy(:,2) - tvy(:,1), tvz(:,2) - tvz(:,1) ];
edge23 = [ tvx(:,3) - tvx(:,2), tvy(:,3) - tvy(:,2), tvz(:,3) - tvz(:,2) ];
normals = cross(edge12, edge23, 2);

% Length of the normals.
nl = sqrt(sum(normals.*normals,2));

% Triangle area is equal to the half of the normal length.
mesh.tri_a = nl/2;

% Now normalize the normals
normals = normals ./ repmat(nl, 1, 3);

mesh.nx = normals(:,1);
mesh.ny = normals(:,2);
mesh.nz = normals(:,3);
