function ndtri = rmduptris(tri)
% ndtri = rmduptris(tri)
%
% Removes duplicated triangles - ones which are referring to the same vertices.
% Removes all the instances of the duplicated triangles, this is different from
% the behavior of the "unique" function which revoves all instances except one.
% The supposed usage is removal of the inner triangles when composing
% a geometry for analysis from simpler primitives.
% 

stri = sort(tri, 2); % sort each row

N = size(tri,1);
[ utri, i, iu ] = unique(stri, 'rows');
nu = accumarray(iu, 1); % number of such triangles for each element of utri
ntri = nu(iu); % number of such triangles for each element of tri
ind = find(ntri == 1); % indices of non-duplicated entries

ndtri = tri(ind,:);
