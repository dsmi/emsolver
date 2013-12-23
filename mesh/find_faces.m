function faces = find_faces(mesh, cx, cy, cz, nx, ny, nz, radius)
% faces = find_faces(mesh, cx, cy, cz, nx, ny, nz, radius)
%
% Finds the mesh faces lying in a particular plane and located within
% a particular radius around a given point. The plane is given by a normal
% and a point lying in the plane, the radius is measured from the same point.
% Typical usage is to find contact faces.
%   Params:
%     mesh       - struct containing the mesh data as returned by init_mesh.
%     cx, cy, cz - point in the plane of interest
%     nx, ny, nz - the plane normal
%     radius     - radius to find edges in
%   Return values:
%     faces      - list of faces found. Row vector.

% Distance from the plane
pr = nx*(mesh.x - cx) + ny*(mesh.y - cy) + nz*(mesh.z - cz);

% Distance from the given point
r = sqrt((mesh.x - cx).^2 + (mesh.y - cy).^2 + (mesh.z - cz).^2);

% Distances for the triangle vertices
tri_pr = pr(mesh.tri);
tri_r  = r(mesh.tri);

% Find the conforming triangles
faces = find(abs(sum(tri_pr,2)) < 1e-8 & abs(tri_r(:,1)) <= radius ...
              & abs(tri_r(:,2)) <= radius & abs(tri_r(:,3)) <= radius).';
