function [tri, x, y, z] = mkdisc(r, n, nr)
% [tri, x, y, z] = mkdisc(r, n, nr)
%
% Makes a flat disc lying in the Y-Z plane. Radius is r, number of edges
% around the disc is n, number of edges along the radius is nr.
% The normal oriented according to the right-hand rule points along
% the positive direction of X-axis.

% Start from creating triangles sharing vertex at the center
zeros_n1 = zeros(1,n+1);
x = zeros_n1;
y = (r/nr)*ones(1,n+1);
z = zeros_n1;

[ x, y, z ] = rotmesh(x, y, z, linspace(0, 2*pi, n+1), zeros_n1, zeros_n1);

% Vertex at the center of the disc
x(n+2) = 0;
y(n+2) = 0;
z(n+2) = 0;

v1 = (1:n);
v2 = v1 + 1;
v3 = repmat(n+2, 1, n);

tri = [ v1', v2', v3' ];

if nr > 1,
    % The rest of the disc is made from a panel coiled around X axis
    [ ptri, px, py, pz ] = mkpanel(2*pi, (r/nr)*(nr-1), n, nr-1);
    [ px, py, pz ] = rotmesh(px, py, pz, 0, pi/2, 0);
    [ px, py, pz ] = move(px, py, pz, 0, (r+r/nr)/2, 0);
    [ px, py, pz ] = rotmesh(px, py, zeros(size(pz)), ...
						   pz+pi, zeros(size(px)), zeros(size(px)));
	
	tri = [ tri; (ptri + length(x)) ];
	x = [ x px ];
	y = [ y py ];
	z = [ z pz ];
end

% Remove duplicated vertices
[ tri, x, y, z ] = rmdups(tri, x, y, z);
