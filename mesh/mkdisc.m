function [tri, x, y, z] = mkdisc(r, n, nr, r0)
% [tri, x, y, z] = mkdisc(r, n, nr, r0)
%
% Makes a flat disc lying in the Y-Z plane. Radius is r, number of edges
% around the disc is n, number of edges along the radius is nr.
% The normal oriented according to the right-hand rule points along
% the positive direction of X-axis.
% r0 is an optional parameter, default is zero, and specifies the radius
% of the hole in the middle. If it is nonzero, the disc becomes a ring
% (or, technically speaking, an annulus)

if ~exist('r0')
    r0 = 0.0;
end

% Start from creating triangles sharing vertex at the center
if 0.0 == r0
    x = zeros(1,n+1);
    y = (r/nr)*ones(1,n+1);
    z = zeros(1,n+1);

    [ x, y, z ] = rotmesh(x, y, z, linspace(0, 2*pi, n+1), z*0.0, z*0.0);

    % Vertex at the center of the disc
    x(n+2) = 0;
    y(n+2) = 0;
    z(n+2) = 0;

    v1 = (1:n);
    v2 = v1 + 1;
    v3 = repmat(n+2, 1, n);

    tri = [ v1', v2', v3' ];

    % Parameters of the rest of the disc
    r1 = (r/nr)*(nr-1);
    nr1 = nr-1;
else
    % No central section
    tri = [ ];
    x = [ ];
    y = [ ];
    z = [ ];

    % Parameters of the rest of the disc
    r1 = r0;
    nr1 = nr;
end

if nr1 > 0,
    % The rest of the disc is made from a panel coiled around X axis
    [ ptri, px, py, pz ] = mkpanel(2*pi, r1, n, nr1);
    [ px, py, pz ] = rotmesh(px, py, pz, 0, pi/2, 0);
    [ px, py, pz ] = move(px, py, pz, 0, (r-r1/2), 0);
    [ px, py, pz ] = rotmesh(px, py, 0.0*pz, pz+pi, 0.0*px, 0.0*px);
	
    tri = [ tri; (ptri + length(x)) ];
    x = [ x px ];
    y = [ y py ];
    z = [ z pz ];
end

% Remove duplicated vertices
[ tri, x, y, z ] = rmdups(tri, x, y, z);
