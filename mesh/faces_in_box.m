function faces = faces_in_box( tri, x, y, z, xmin, ymin, zmin, xmax, ymax, zmax)
% faces = faces_in_box( tri, x, y, z, xmin, ymin, zmin, xmax, ymax, zmax )
%
% Finds the mesh faces lying inside the given box. All edges of the face
% needs to be in the box for the box to be considered in.
% Lower limit is included, upper is not.
%   Params:
%     tri, x, y, z                        - triangles to examine
%     xmin, ymin, zmin, xmax, ymax, zmax  - box
%   Return values:
%     faces      - list of faces found. Row vector.

x123 = x( tri );
y123 = y( tri );
z123 = z( tri );

xin = all( ( x123 >= xmin ) & ( x123 < xmax ), 2 );
yin = all( ( y123 >= ymin ) & ( y123 < ymax ), 2 );
zin = all( ( z123 >= zmin ) & ( z123 < zmax ), 2 );

faces = transpose( find( xin & yin & zin ) );
