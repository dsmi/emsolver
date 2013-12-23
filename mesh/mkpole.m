function [tri, x, y, z] = mkpole(l, r, nl, n, nr)
% [tri, x, y, z] = mkpole(l, r, nl, n, nr)
%
% Makes a pole - hollow cylinder parallel to X axis with discs at both ends.
% The resulting triangles are oriented CCW when looking from outside,
% so the normal oriented according to the right-hand rule is pointing
% outwards.
%  l  - length of the pole
%  r  - radius of the pole
%  nl - number of edges along the pole
%  n  - number of edges around the cross section
%  nr - number of edges along the radius in the discs at the ends

% Start from a panel perpendicular to Z
[ tri, x, y, z ] = mkpanel(l, 2*pi, nl, n);

% Bend this panel around x. Minus is to get the normals looking outside.
[ x, y, z ] = rotate(x, r*ones(size(y)), zeros(size(z)), ...
                               -y+pi, zeros(size(y)), zeros(size(y)));

% Add discs at the ends
[ dtri, dx, dy, dz ] = mkdisc(r, n, nr);
[ dx, dy, dz ] = move(dx, dy, dz, l/2, 0, 0);
	
tri = [ tri; (dtri + length(x)) ];
x = [ x dx ];
y = [ y dy ];
z = [ z dz ];

[ dtri, dx, dy, dz ] = mkdisc(r, n, nr);
%dtri = [ dt(:,2) dt(:,1) dt(:,3) ]; % To get the correct orientation.
[ dx, dy, dz ] = rotate(dx, dy, dz, 0, pi, 0);
[ dx, dy, dz ] = move(dx, dy, dz, -l/2, 0, 0);
	
tri = [ tri; (dtri + length(x)) ];
x = [ x dx ];
y = [ y dy ];
z = [ z dz ];


% Remove duplicated vertices
[ tri, x, y, z ] = rmdups(tri, x, y, z);

