
addpath(genpath([ pwd, '/..' ]));

mesh_globals;

% The source primitive - two loops, different sizes and resolutions
[ tri1, x1, y1, z1 ] = mkring(1, 0.1, 3, 3);
[ tri2, x2, y2, z2 ] = mkring(0.1, 0.01, 3, 3);
m_x = [ x1 x2 ];
m_y = [ y1 y2 ];
m_z = [ z1 z2 ];
m_tri = [ tri1; (tri2 + length(x1)) ];

init_edges;
init_mesh;

% Find edges which belong to each primitive
trie1 = m_tri_edges(1:size(tri1)(1),:);
edges1 = unique(trie1(:));
trie2 = m_tri_edges(size(tri1)(1)+1:end,:);
edges2 = unique(trie2(:));
m_objedges = { edges1 edges2 };

%[ m_tri, m_x, m_y, m_z ] = mkring(1, 0.1, 3, 3);

% Number of edges
ne = size(m_edges,1)


[ IL, IS ] = mkloopstar;
size_IL = size(IL)
size_IS = size(IS)
T = [ IL IS ];
size_T = size(T)
rank_T = rank(T)

% The arbitrary loop basis function coefficients.
X = (1:size(IL, 2))';

X = IL*X;
D = mkdivmat;
face_div = D*X;

max(face_div)
