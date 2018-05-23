function [ IL, IS ] = mkloopstar(mesh, edges)
% [ IL, IS ] = mkloopstar(mesh, edges)
%
% Makes IL and IS matrices which transform RWG basis functions to loops and
% stars. The matrices are defined as:
%    IL*Il + IS*Is = I, [L] = IL'*[f], [S] = IS'*[f]
% where Il and Is are the loops and stars basis function coefficients
% correspondingly, [f] is conventional RWG basis, [L] is the loop basis
% and [S] is the star basis. The matrix dimensions are IL(NxNL) and IS(NxNS),
% NL is number of loops and NS is number of stars.
% There is an exapmle how EFIE for PEC looks like when using loop-star basis:
% T = [ IL IS ]     - transform matrix
% Z*I=V
%   Zls = Tt*Z*T  (Z = <<G,f_n>,f_m>, f_n - expansion, f_m - testing)
%   I = T*Ils     (I is the coefficients vector of basis f_n)
%   Vls = Tt*V    (V = <Einc,f_m>, excitation)
%

% If the edges is not defined then create loops and stars
% for the whole mesh, separately for each of the shapes.
if ~exist('edges')
    % for each shape separately, and glue together
    [ IL, IS ] = mkloopstar(mesh, mesh.shape_edges{1});
    for i=2:length(mesh.shape_edges)
	[ il, is ] = mkloopstar(mesh, mesh.shape_edges{i});
	IL = [ IL il ];
	IS = [ IS is ];
    end
    
    return
end

% Triangles of an edge, global triangle indices
edge_tris = mesh.edge_tris(edges,:); 

% triangles of the current shape
tris = unique(edge_tris(:));

% global-to-local triangle index
tris_g2l = zeros(size(mesh.tri, 1), 1);
tris_g2l(tris) = 1:length(tris);

% Triangles of an edge, local triangle indices
edge_tris_local = tris_g2l(edge_tris); 

% edges for a triangle and signs
tri_edges = mesh.tri_edges(tris,:);
tri_edges_s = mesh.tri_edges_s(tris,:);

% vertices of an edge, global vertex indices
edge_verts = mesh.edges(edges,:); 

% vertices of the current shape
verts = unique(edge_verts(:));

% global-to-local vertex index
verts_g2l = zeros(length(mesh.x), 1);
verts_g2l(verts) = 1:length(verts);

% vertices of an edge, local vertex indices
edge_verts_local = verts_g2l(edge_verts); 

% edge lengths
edge_l = mesh.edge_l(edges,:);

% Number of edges in the shape
ne = size(edge_tris, 1);

% Number of triangles in the shape
nt = size(tris, 1);

% Number of vertices in the shape
nv = length(verts);

% Number of rows in loop and star matrices - total number of edges
nrows = size(mesh.edges, 1);
	
% Number of loops
nl = nv - 1;

% Loops. We are going to create as many loops as the number of vertices
% and delete extra loops later.
% An edge has two adjoining faces, positive one which has the same
% traversal direction as the edge, and negative one with the opposite
% direction. The basis function is directed from positive triangle to
% negative one. Let the loops be oriented CCW when looking from outside,
% like the faces are, then for the loop anchored to first vertex this edge
% has positive direction and negative for the second vertex.
m = repmat(edges, 1, 2);
n = edge_verts_local;
val = repmat([ 1 -1 ], ne, 1)./repmat(edge_l,1,2);
IL = zeros(nrows, nv);
IL(sub2ind(size(IL), m, n)) = val;
%IL = sparse(m, n, val, nrows, nv);

% Remove extra loops.
IL(:,nl+1:end) = [];

if rank(IL) ~= size(IL,2)
    error('mkloopstar: vertex-anchored loops are incorrect - check mesh.shape_edges');
end

% The loop/cycle basis created above is complete if the shape is simple 
% (has no holes), if it is not then a number of additional base cycles
% needs to be found. Determine how many cycles we need to create based
% on the Euler characteristic of the shape.
naddl = 2 - (nv - ne + nt);
if naddl > 0

    nl = nl + naddl;

    % de Pina algorithm is used to find the additional cycles.

    % spanning tree and vector representation
    [ se, sei ] = mst(edge_tris_local);
    vse = zeros(length(se), ne);
    vse(sub2ind(size(vse), 1:length(se), sei')) = 1;
    
    while size(IL, 2) < nl

	% vector representation of the cycles
	%size(IL)
	cycles = +full(IL(edges,:) ~= 0);

	% compute vector orthogonal to all the cycles
	% and to the edges of the spanning tree
	[ r2, orth ] = modrank2( [ cycles' ; vse ] );

	% next, we need to find the cheapest cycle with <orth,C> ~= 0
	% to do that, build a graph with two copies of nodes and
	% edges, n+ n- and e+ e-. edges from orth change sides, others don't
	pose = edge_tris_local      + [ orth*0, orth*nt ];
	nege = edge_tris_local + nt - [ orth*0, orth*nt ];
	pngraph = [ pose ; nege ];
	
	% for all n find shortest path from n+ to n-
	mind = 1e99;
	for i=1:nt
	    [ d, p ] = bfs(pngraph, i, i+nt);
	    if d(i+nt) < mind
		mind = d(i+nt);
		pred = p;
		v = i+nt;
	    end
	end

	% build the shortest path just found
	path = [];
	while v ~= pred(v)
	    path = [ path; [ v, pred(v) ] ];
	    v = pred(v);
	end

	% then, based on the path, determine the cycle edges
	path_edges_v = ismember(sort(pngraph, 2), sort(path, 2), 'rows');
	cycle_edges_v = mod(path_edges_v(1:ne)+path_edges_v(ne+1:end), 2);

	% indices of the edges forming the cycle
	cycle_edges = find(cycle_edges_v);
	
	% finally, we need to determine the orientation of the edges
	% in the cycle.
	cycle_edges_s = ones(length(cycle_edges), 1);
	cet = edge_tris_local(cycle_edges,:);
	v = cet(1,2);
	cet(1,:) = [ 0 0 ];
	while any(cet)
	    idx = find(v == cet(:,1));
	    if ~isempty(idx)
		v = cet(idx,2);
	    else
		idx = find(v == cet(:,2));
		cycle_edges_s(idx) = -1;
		v = cet(idx,1);
	    end
	    cet(idx,:) = [ 0 0 ];
	end
	
	% And add the cycle to the matrix!
	m = edges(cycle_edges);
	n = ones(size(cycle_edges));
	val = cycle_edges_s./edge_l(cycle_edges);
	il = zeros(nrows, 1);
	il(sub2ind(size(il), m, n)) = val;
	% il = sparse(m, n, val, nrows, 1);
	IL = [ IL, il ];

    end
end

% Number of stars
ns = ne - nl;

% Stars. We are going to create as many stars as the number of triangles
% and delete extra stars later.
m = tri_edges;
n = repmat((1:nt)', 1, 3);
pn = [ 1 -1 ];
val = pn(tri_edges_s)./mesh.edge_l(tri_edges);
IS = zeros(nrows, nt);
IS(sub2ind(size(IS), m, n)) = val;
%% IS = sparse(m, n, val, nrows, nt);

% Remove extra stars.
IS(:,ns+1:end) = [];
