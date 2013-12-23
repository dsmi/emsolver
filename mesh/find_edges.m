function [edges, edges_s] = find_edges(mesh, cx, cy, cz, nx, ny, nz, radius)
% [edges, edges_s] = find_edges(mesh, cx, cy, cz, nx, ny, nz, radius)
%
% Finds the mesh edges lying in a particular plane and located within
% a particular radius around a given point. The plane is given by a normal
% and a point lying in the plane, the radius is measured from the same point.
% Typical usage is to find excitation edges.
%   Params:
%     mesh       - struct containing the mesh data as returned by init_mesh.
%     cx, cy, cz - point in the plane of interest
%     nx, ny, nz - the plane normal
%     radius     - radius to find edges in
%   Return values:
%     edges      - list of edges found. Row vector.
%     edges_s    - edge signs, 1 or -1; vector of the same size as edges.
%                  Becuase the basis function is directed from positive
%                  triangle to negative, the edge is positive if the vector
%                  from positive triangle free vertex to negative one is
%                  directed along the given normal (dot product is positive).

% Distance from the plane
pr = nx*(mesh.x - cx) + ny*(mesh.y - cy) + nz*(mesh.z - cz);

% Distance from the given point
r = sqrt((mesh.x - cx).^2 + (mesh.y - cy).^2 + (mesh.z - cz).^2);

% Distances for the edge points.
pr1 = pr(mesh.edges(:,1));
pr2 = pr(mesh.edges(:,2));
r1  = r(mesh.edges(:,1));
r2  = r(mesh.edges(:,2));

% Find the conforming edges.
edges = find(abs(pr1+pr2) < 1e-8 & abs(r1) <= radius & abs(r2) <= radius);

% Vector from positive triangle free vertex to negative for the edges found
pos_to_neg_x = mesh.x(mesh.free_vert(edges, 2)) - mesh.x(mesh.free_vert(edges, 1));
pos_to_neg_y = mesh.y(mesh.free_vert(edges, 2)) - mesh.y(mesh.free_vert(edges, 1));
pos_to_neg_z = mesh.z(mesh.free_vert(edges, 2)) - mesh.z(mesh.free_vert(edges, 1));

% Determine the sign by calculating the dot product
edges_s = sign(pos_to_neg_x.*nx + pos_to_neg_y.*ny + pos_to_neg_z.*nz);
