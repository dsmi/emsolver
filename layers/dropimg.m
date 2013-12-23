function [ a, b ] = dropimg(a0, b0)
% [ a, b ] = dropimg(a0, b0)
%
% Drops images with negligible contribution. In the current implementation
% drops the images which are ten orders of magnitude less than the maximum
% one.
% TODO: think of better way to filter the images.
%

a = a0;
b = b0;

if ~isempty(a0)
	
	% magnitudes of the images
	magn = abs(a0.*exp(-j*sqrt(b0.*b0)));
	maxm = max(magn);
	to_drop = find(magn < maxm*1e-10);
	
	a(to_drop) = [];
	b(to_drop) = [];

end
