function [ e, v ] = mkcir2d(r, n)
% [ e, v ] = mkcir2d(r, n)
%
% Makes a circle centered at the origin. Radius is r, number of edges
% around the circle is n. The resulting circle is oriented CCW as it
% is required for outer boundary, to change the orientantion just
% scale it with either of the coefficients set to -1.
%

a = linspace(0, 2*pi, n+1).';
vx = cos(a)*r;
vy = sin(a)*r;
v = [ vx vy ];

e = [ (1:n)' (2:n+1)' ];

% Remove duplicated vertices
[ e, v ] = rmdups2d(e, v);
