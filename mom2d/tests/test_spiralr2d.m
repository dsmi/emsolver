function test_spiralr2d
%
% Compute resistance of a spiral.
%

% The geometry - a spiral
nt = 3;
w = 0.1;
maxl = 0.1;
[ e, v, l ] = mkspiral2d(nt,w,maxl);

% Find ports
p1 = find_edges2d(e, v, -0.5, 0, w/2*1.1);
p2 = find_edges2d(e, v, -2, -1.5, w/2*1.1);
ports = { p1' p2' };

%plotmesh2d(e,v,ports,0);

Y = extracty2(e, v, ports, @intg_lapsl2d, @intg_lapdl2d);

R = Y(1,1);
R_test = w/l;

assertEquals(R_test, R, R_test*5e-2);
