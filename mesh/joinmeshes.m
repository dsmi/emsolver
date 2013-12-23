function [ tri, x, y, z ] = joinmeshes(ctri, cx, cy, cz)
% [ tri, x, y, z ] = joinmeshes(ctri, cx, cy, cz)
%
% Joins together a few meshes, updates the vertex indices of the second and
% the subsequent meshes correspondingly. ctri, cx, cy, cz are the cell arrays
% of the same length.
%


x = cell2mat(reshape(cx, 1, length(cx)));
y = cell2mat(reshape(cy, 1, length(cy)));
z = cell2mat(reshape(cz, 1, length(cz)));

% index offsets
o = cumsum(cellfun('length', cx));
o = num2cell([ 0 o(1:end-1) ]);

t = cellfun(@(t,o) t+o, ctri, o, 'UniformOutput', false);

tri = cell2mat(reshape(t, length(t), 1));
