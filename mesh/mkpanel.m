function [tri, x, y, z] = mkpanel(xlen, ylen, xnum, ynum)
% [tri, x, y, z] = mkpanel(xlen, ylen, xnum, ynum)
%
% Makes a flat rectangular panel lying in the X-Y plane. Dimensions of
% the panel are (-xlen/2 - xlen/2), (-ylen/2 - ylen/2) along x and y
% correspondingly. The panel consists of rectangular segments each
% consisting of two triangles, the number of segments are xnum and ynum
% along x and y correspondingly. The resulting triangles are oriented CCW,
% so the normal oriented according to the right-hand rule points along
% the positive direction of Z-axis.

xgrid = linspace(-xlen/2, xlen/2, xnum + 1);
x = xgrid(ones(1,ynum + 1),:);
x = x(:)';

ygrid = linspace(-ylen/2, ylen/2, ynum + 1)'; % Notice the transpose.
y = ygrid(:,ones(xnum + 1,1));
y = y(:)';

z = zeros(1,length(x));

% Number of rectangular segments.
num_of_seg = xnum*ynum;

% Vertices of the segment
seg_ind = (1:num_of_seg);
v1 = seg_ind + floor((seg_ind - 1)./ynum);
v4 = v1 + 1;
v2 = v1 + ynum + 1;
v3 = v2 + 1;

% Number of triangles, each segment consists of the two
num_of_tri = num_of_seg * 2;

tri = zeros(num_of_tri,3);
tri(1:2:num_of_tri,:) = [ v1', v2', v3' ];
tri(2:2:num_of_tri,:) = [ v3', v4', v1' ];
