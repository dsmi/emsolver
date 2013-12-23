function test_stripr2d
%
% Compute resistance of a strip.
%

% The geometry - a strip
length = 100;
width = 1;
nl = 200;
nw = 2;
helm_k = 1e-10;
[ e, v ] = mkrect2d(length, width, nl, nw);

% Find ports
p1 = find_edges2d(e, v, -length/2, 0, width/2*(1+1e-10));
p2 = find_edges2d(e, v, length/2, 0, width/2*(1+1e-10));
ports = { p1' p2' };

Y = extracty2(e, v, ports, @intg_lapsl2d, @intg_lapdl2d);

R = Y(1,1);
R_test = width/length;

assertEquals(R_test, R, R_test*1e-2);
