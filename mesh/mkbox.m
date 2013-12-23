function [tri, x, y, z, c1, c2] = mkbox(sx, sy, sz, nx, ny, nz, panels)
% [tri, x, y, z, c1, c2] = mkbox(sx, sy, sz, nx, ny, nz, panels)
%
% Makes a 3D box. Dimensions of the box are (-sx/2 - sx/2), (-sy/2 - sy/2),
% (-sz/2 - sz/2) along x, y and z correspondingly; mesh resolution is nx,
% ny, nz.
% The optional parameter panels allow to skip certain of the panels
% forming the box, thus creating a tube of rectangular crossection, or a
% corner and so on. If the panels parameter is defined, it is treated as
% a sum of flag values, which are defined as follows:
%   1  - panel perpendicular to X axis, located at -sx/2
%   2  - panel perpendicular to X axis, located at sx/2
%   4  - panel perpendicular to Y axis, located at -sy/2
%   8  - panel perpendicular to Y axis, located at sy/2
%   16 - panel perpendicular to Z axis, located at -sz/2
%   32 - panel perpendicular to Z axis, located at sz/2
% If a flag is not set, the corresponding panel is ommited.
% The resulting triangles are oriented CCW when looking from outside,
% so the normal oriented according to the right-hand is pointing outwards.
% Outputs c1 and c2 are the lists of indices of triangles forming edges of
% the box which are perpendicular to x axis.

% First panel perpendicular to Z axis
[ pz1_tri, pz1_x, pz1_y, pz1_z ] = mkpanel(sx, sy, nx, ny);
[ pz1_x, pz1_y, pz1_z ] = move(pz1_x, pz1_y, pz1_z, 0, 0, sz/2);

% Second panel perpendicular to Z axis
[ pz2_tri, pz2_x, pz2_y, pz2_z ] = mkpanel(sx, sy, nx, ny);
[ pz2_x, pz2_y, pz2_z ] = rotate(pz2_x, pz2_y, pz2_z, 0, pi, 0);
[ pz2_x, pz2_y, pz2_z ] = move(pz2_x, pz2_y, pz2_z, 0, 0, -sz/2);


% First panel perpendicular to X axis
[ px1_tri, px1_x, px1_y, px1_z ] = mkpanel(sz, sy, nz, ny);
[ px1_x, px1_y, px1_z ] = rotate(px1_x, px1_y, px1_z, 0, pi/2, 0);
[ px1_x, px1_y, px1_z ] = move(px1_x, px1_y, px1_z, sx/2, 0, 0);

% Second panel perpendicular to X axis
[ px2_tri, px2_x, px2_y, px2_z ] = mkpanel(sz, sy, nz, ny);
[ px2_x, px2_y, px2_z ] = rotate(px2_x, px2_y, px2_z, 0, -pi/2, 0);
[ px2_x, px2_y, px2_z ] = move(px2_x, px2_y, px2_z, -sx/2, 0, 0);


% First panel perpendicular to Y axis
[ py1_tri, py1_x, py1_y, py1_z ] = mkpanel(sx, sz, nx, nz);
[ py1_x, py1_y, py1_z ] = rotate(py1_x, py1_y, py1_z, -pi/2, 0, 0);
[ py1_x, py1_y, py1_z ] = move(py1_x, py1_y, py1_z, 0, sy/2, 0);

% Second panel perpendicular to Y axis
[ py2_tri, py2_x, py2_y, py2_z ] = mkpanel(sx, sz, nx, nz);
[ py2_x, py2_y, py2_z ] = rotate(py2_x, py2_y, py2_z, pi/2, 0, 0);
[ py2_x, py2_y, py2_z ] = move(py2_x, py2_y, py2_z, 0, -sy/2, 0);

% Default value of the panels optional parameter
if ~exist('panels'),
	panels = 1+2+4+8+16+32;
endif

% Merge all together
tri = [];
x = [];
y = [];
z = [];

if rem(panels,64) >= 32, % Testing for the 32 flag
	tri = [ tri; (pz1_tri + length(x)) ];
	x = [ x pz1_x ];
	y = [ y pz1_y ];
	z = [ z pz1_z ];
endif

if rem(panels,32) >= 16, % Testing for the 16 flag
	tri = [ tri; (pz2_tri + length(x)) ];
	x = [ x pz2_x ];
	y = [ y pz2_y ];
	z = [ z pz2_z ];
endif

if rem(panels,4) >= 2, % Testing for the 2 flag
	tri = [ tri; (px1_tri + length(x)) ];
	x = [ x px1_x ];
	y = [ y px1_y ];
	z = [ z px1_z ];
endif

if rem(panels,2) >= 1, % Testing for the 1 flag
	tri = [ tri; (px2_tri + length(x)) ];
	x = [ x px2_x ];
	y = [ y px2_y ];
	z = [ z px2_z ];
endif

if rem(panels,16) >= 8, % Testing for the 8 flag
	tri = [ tri; (py1_tri + length(x)) ];
	x = [ x py1_x ];
	y = [ y py1_y ];
	z = [ z py1_z ];
endif

if rem(panels,8) >= 4, % Testing for the 4 flag
	tri = [ tri; (py2_tri + length(x)) ];
	x = [ x py2_x ];
	y = [ y py2_y ];
	z = [ z py2_z ];
endif

% Remove duplicated vertices
[ tri, x, y, z ] = rmdups(tri, x, y, z);
