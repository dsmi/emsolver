function [ tri x y z ] = mkwireddisk( t, r, l1, l2, rw, nl, n, nr, nw, nrw )
% [ tri x y z ] = mkwireddisk( t, r, l1, l2, rw, nl, n, nr, nw, nrw )
%
% Makes a thick disk (cylinter) with two wires along the axis of symmetry.
%
% The resulting triangles are oriented CCW when looking from outside,
% so the normal oriented according to the right-hand rule is pointing
% outwards.
%
%  t  - trickness of the disk
%  r  - radius of the disk
%  l1, l2 - length of the upper (x < 0) and lower wires correspondingly
%  rw - radius of the wires
%  nl - number of edges along the disk edge
%  n  - number of edges around the cross section
%  nr - number of edges along the radius of the top/bottom
%  nw - number of edges along the wire
%  nrw - number of edges along the radius of the wire
%

%% [ tri, x, y, z ] = mkpole(t, r, nl, n, nr, rw);

% Outer surface of the disk
[ tri, x, y, z ] = mktube(t, r, nl, n);

% Add top and bottom covers
[ dtri, dx, dy, dz ] = mkdisc(r, n, nr, rw);
[ dx, dy, dz ] = move(dx, dy, dz, t/2, 0, 0);

[ tri, x, y, z ] = joinmeshes( { tri, dtri }, { x, dx }, { y, dy }, { z, dz } );

[ dtri, dx, dy, dz ] = mkdisc(r, n, nr, rw);
[ dx, dy, dz ] = rotmesh(dx, dy, dz, 0, pi, 0); % to have outward normals
[ dx, dy, dz ] = move(dx, dy, dz, -t/2, 0, 0);
	
[ tri, x, y, z ] = joinmeshes( { tri, dtri }, { x, dx }, { y, dy }, { z, dz } );

% Top wire (x < 0)
if l1 ~= 0
    [ wt, wx, wy, wz ] = mktube( l1, rw, nw, n );
    [ wx, wy, wz ] = move(wx, wy, wz, -t/2 - l1/2, 0, 0);

    [ tri, x, y, z ] = joinmeshes( { tri, wt }, { x, wx }, { y, wy }, { z, wz } );
end

% Bottom wire (x > 0)
if l2 ~= 0
    [ wt, wx, wy, wz ] = mktube( l2, rw, nw, n );
    [ wx, wy, wz ] = move(wx, wy, wz, t/2 + l2/2, 0, 0);

    [ tri, x, y, z ] = joinmeshes( { tri, wt }, { x, wx }, { y, wy }, { z, wz } );
end

% Covers of the wires
[ dtri, dx, dy, dz ] = mkdisc(rw, n, nrw);
[ dx, dy, dz ] = move(dx, dy, dz, t/2 + l2, 0, 0);

[ tri, x, y, z ] = joinmeshes( { tri, dtri }, { x, dx }, { y, dy }, { z, dz } );

[ dtri, dx, dy, dz ] = mkdisc(rw, n, nrw);
[ dx, dy, dz ] = rotmesh(dx, dy, dz, 0, pi, 0); % to have outward normals
[ dx, dy, dz ] = move(dx, dy, dz, -t/2 - l1, 0, 0);
	
[ tri, x, y, z ] = joinmeshes( { tri, dtri }, { x, dx }, { y, dy }, { z, dz } );
    
% Remove duplicated vertices
[ tri, x, y, z ] = rmdups(tri, x, y, z);
