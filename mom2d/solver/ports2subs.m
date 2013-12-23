function [ faceidx portidx ] = ports2subs(ports)
% [ faceidx portidx ] = ports2subs(ports)
%
% Given a cell vector of vectors, each of which lists faces (or edges)
% forming a port, this function builds two vectors of indices of the same
% size. First one is a concatenation of all vectors of faces, and the
% second vector with indices of the corresponding ports.
% For example, if the following cell array is given:
%   ports = { [ 2 4 ] [ 7 8 9 ] }
% The results are:
%   faceidx = [ 2 4 7 8 9 ]
%   portidx = [ 1 1 2 2 2 ]
%
% Inputs:
%   ports   - cell vector of vectors, each of which lists faces belonging
%             to a port. May be whether row or column vector, but if it is
%             a row vector the nested vectors must be row ones as well, and
%             if it is column the nested ones must be columns.
% Outputs:
%   faceidx - row vector of all faces belonging to ports
%   portidx - row vector of the ports indices of each face
%

faceidx = cell2mat(ports);
faceidx = reshape(faceidx, 1, length(faceidx));

p = cellfun('length', ports);
% Here is what the trick below does - given vector p, it builds vector
% portidx of length sum(p) such that the first p(1) elements are 1, the
% next p(2) elements are 2 and so on.
fnz = find(p > 0);
a = zeros(max(p),length(p));
a(sub2ind(size(a),p(fnz),fnz)) = fnz;
b = flipud(cumsum(flipud(a),1));
portidx = b(b > 0);
portidx = reshape(portidx, 1, length(portidx));
