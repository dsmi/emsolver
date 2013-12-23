function [ e, v, l ] = mkspiral2d(nt,w,maxl)
% [ e, v, l ] = mkspiral2d(nt,w)
%
% Inputs:
%   nt   - number of turns
%   w    - width
%   maxl - maximum allowed edge length. Too long edges are
%          splitted.
% Outputs:
%   e, v - edges and vertices
%   l    - length of the spiral

% Length of each particular step
sl = reshape([ (1:nt)-0.5 ; 1:nt ; 1:nt ; (1:nt)+0.5 ], nt*4, 1);

x0 = -0.5;
y0 = 0;
dx = repmat([ 0 1 0 -1 ]', nt, 1);
dy = repmat([ 1 0 -1 0 ]', nt, 1);

% Centers of the bends
xc = [ x0 ; x0 + cumsum(sl.*dx) ];
yc = [ y0 ; y0 + cumsum(sl.*dy) ];
l = sum(abs(sl.*dx)) + sum(abs(sl.*dy));

% Create bends
dxb = [ repmat([ -w/2 -w/2 w/2 w/2 ]', nt, 1); 0 ];
dyb = [ repmat([ -w/2 w/2 w/2 -w/2 ]', nt, 1); -w/2 ];
dyb(1) = 0;
x1 = xc + dxb;
y1 = yc + dyb;
x2 = xc - dxb;
y2 = yc - dyb;

v = [ [ x1; x2 ]  [ y1; y2 ] ];

nb = nt*4; % Number of bends
e = [ [ 1 nb+2 ];     [ (nb+2:nb*2+1)' (nb+3:nb*2+2)' ]; ...
      [ nb*2+2 nb+1]; [ (nb+1:-1:2)'   (nb:-1:1)'     ] ];

% Edge beginnings and ends
r1 = v(e(:,1), :);
r2 = v(e(:,2), :);

% Calculate edge lengths.
el = sqrt(sum((r2 - r1).^2,2)); % Column vector of length N.

% Split the edges if too long. Vectorize this?
newe = [];
for ie=1:size(e,1),
	n = ceil(el(ie)/maxl-1e-6);
	if n > 1,
		x = linspace(r1(ie,1), r2(ie,1), n+1);
		y = linspace(r1(ie,2), r2(ie,2), n+1);
		newe = [ newe; size(v,1) + [ (1:n)' (2:n+1)' ] ];
		v = [ v; [ x' y' ] ];
	else
		newe = [ newe; e(ie,:) ];
	endif
end

e = newe;

% Remove duplicated vertices, which may be created when splitting edges
[ e, v ] = rmdups2d(e, v);
