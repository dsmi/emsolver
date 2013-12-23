function [ tri_e_over_r, edge_f_e_over_r, edge_e_over_r, ...
            edge_f_cross_grad ] = calc_edge_intg(mesh, k, edges, chunk_size)
% [ tri_e_over_r, edge_f_e_over_r, edge_e_over_r, ...
%           edge_f_cross_grad ] = calc_edge_intg(mesh, k, edges, chunk_size)
%
% !!! This function is now deprecated and is used for testing only.
%
% Computes potential integrals associated with pairs of edges.
% The integrals evaluator (see i_calc) is suited for computing the integrals
% for face-observation point pairs, while indeed we are interested in the
% integrals associated with pairs of edges. A pair of edges consists of
% a testing (or observation) edge and a source edge. Each integral associated
% with a pair of edges invloves integration over two faces: positive and
% negative ones of the source edge and two evaluation points which are
% centers of the positive and negative faces of the testing edge.
% This function also applies the basis function normalization multipliers
% which are omitted by the integrals evaluator.
% Given the fully initialized mesh this function prepares input
% for the integrals evaluator, runs it and based on its results fills
% the output arrays.
%
% mesh - struct containing the mesh data as returned by init_mesh.
%
% k - Wavenumber
%
% edges - Allows to specify a subset of edges to evaluate the integrals
%     for. A column array of the edge indices. If integrals for all edges
%     are needed, this vector should be 1:num_of_edges (or empty/omitted).
%     Otherwise integrals are evaluated for pairs formed from edges from
%     this vector. In this case, num_of_edges in the return values below
%     is a number of edges in this vector.
%
% chunk_size - The integrals are computed using the integrals evaluator
%     (see i_calc) which can not handle all the face pairs at once due
%     to the memory limitations if the mesh is big enough, so the entire
%	  face/observation pairs array is being passed to the evaluator by
%     chunks. This variable defines the chunk size. If not defined,
%     a default value of 20K is used, which works ok on my 2GB laptop.
%
% tri_e_over_r - Integrals of the exp(-j*k*R)/R for all the triangle pairs.
%     num_of_tris-by-num_of_tris array. The indices are the testing
%     and source triangles correspondingly, integrals are evaluated at the
%     centers of the testing triangles.
%
% edge_f_e_over_r - Integrals of the f*exp(-j*k*R)/R for all the edge pairs;
%     see the note below for some additional details.
%     num_of_edges-by-num_of_edges-by-2-by-3 array. First two indices are the
%     testing and source edges correspondingly, second index is
%     positive/negative faces of the tesing edge, fourth index is X-Y-Z.
%       The following applies to edge_f_e_over_r, edge_e_over_r and
%     edge_f_cross_grad: each entry of the array corresponds to a pair of
%     edges: edge m which is testing (or observation) edge and edge n
%     which is source edge. The integration is performed over the positive
%     and negative faces of the source edge, the evaluation points are
%     centers of the positive and negative faces of the testing edge.
%     Unlike the integrals evaluator (see i_calc) where the basis function
%     f is not normalized here f is defined as:
%        f = (ln/2An)*(r-r(l))  if r is in positive face
%        f = (ln/2An)*(r(l)-r)  if r is in negative face
%     ln is the length of the edge n and An is the area of the corresponding
%     positive or negative face.
%   
% edge_e_over_r - Integrals of the div(f)*exp(-j*k*R)/R for all the edge
%     pairs. The div(f) is constant in each face and is given by:
%        div(f) = ln/An    if r is in positive face
%        div(f) = -ln/An   if r is in negative face
%     num_of_edges-by-num_of_edges-by-2 array. First two indices are the
%     testing and source edges correspondingly, second index is
%     positive/negative faces of the tesing edge.
%
% edge_f_cross_grad - Integrals of the f X grad(exp(-j*k*R)/R) for all the
%     edge pairs. Notice that R here is the distance from integration point
%     to the observation point, R=R(r,r')=|r-r'|, where r is the obseravtion
%     coordinate and r' is the integration coordinate, and the derivative
%     is taken with respect to the interation coordinate r'.
%     num_of_edges-by-num_of_edges-by-2-by-3 array. First two indices are the
%     testing and source edges correspondingly, second index is
%     positive/negative faces of the tesing edge, fourth index is X-Y-Z.

% If edges is defined then calculate integrals for the given edges only
if exist('edges') && ~isempty(edges),
	nedges = length(edges); % Number of the edges
	% Indices of the triangles used by the edges we need to consider.
    etris = mesh.edge_tris(edges,:);
	trii = unique(etris(:));
	tri = mesh.tri(trii,:);
	% Given the triangle indes in the full triangles vector newtrii
	% vector allows to find index in the subset of triangles used
	% by edges of interest.
	newtrii = zeros(size(mesh.tri,1),1);
	newtrii(trii) = (1:length(trii))';
	edge_tris = newtrii(mesh.edge_tris(edges,:));
	edge_l = mesh.edge_l(edges);
	tri_a = mesh.tri_a(trii);
	free_vert_loc = mesh.free_vert_loc(edges,:);
else
	% Use all otherwise.
	tri = mesh.tri;
	trii = (1:size(tri,1))';
	% Number of the edges
	nedges = size(mesh.edges,1);
	edge_tris = mesh.edge_tris;
	edge_l = mesh.edge_l;
	tri_a = mesh.tri_a;
	free_vert_loc = mesh.free_vert_loc;
end

% Here we start preparing the integrals evaluator input data.

% Number of triangles involves
ntri = size(tri,1);

% The source triangles, ntri-by-3-by-3 array, three points
% per each triangle, the last index is X-Y-Z
tri_r = cat(3, mesh.x(tri), mesh.y(tri), mesh.z(tri));

% Observation points - the triangle centers, ntri-by-3 array,
% the second index is X-Y-Z.
o_r = [ mesh.cx(trii) mesh.cy(trii) mesh.cz(trii) ];

% Now build all the source-observation triangle combinations. The result is
% the matrix with row number corresponding to the observation point and column
% number corresponding to the source triangle.
tri_r = repmat(permute(tri_r, [ 4 1 2 3 ]), ntri, 1);
o_r = repmat(permute(o_r, [ 1 3 2 ]), 1, ntri);

% The integrals evaluator expects vectors of sources/observations,
% so we linearize this matrix column-wise.
% Given the observation point index m and the source point index n one can
% obtain the corresponding observation-source index in the linearized vectors
% as: sub2ind([ntri ntri],m,n).
i_r = reshape(tri_r, [], 3, 3);
i_obs_r = reshape(o_r, [], 3);

% Integrals evaluator uses lots of intermediate data and can not handle
% all face-observation pairs at once due to the memory limitations, so
% we evaluate the integrals by chunks.

% Total number of face/observation pairs
npairs = size(i_r, 1);

f_e_over_r = zeros(3,npairs,3);
e_over_r = zeros(1,npairs);
f_cross_grad = zeros(3,npairs,3);

% If i_chunk_size is defined use it, otherwise use default value
if exist('chunk_size') && ~isempty(chunk_size) && chunk_size > 0,
	chunk_size_to_use = chunk_size;
else
	chunk_size_to_use = 20000;
end

if npairs > 20000,
	fprintf(1, 'Evaluating integrals');
end

last_pair = 0;
while last_pair < npairs,

	first_pair = last_pair + 1;
	last_pair = min(last_pair + chunk_size_to_use, npairs);
	
	pairs_range = first_pair:last_pair;
	ir = i_r(pairs_range,:,:);
	ior = i_obs_r(pairs_range,:);

	% Run the integrals evaluator
	[ fer, er, fxg ] = i_calc(k, ir, ior, 32);
	
	f_e_over_r(:,pairs_range,:)   = fer;
	e_over_r(pairs_range)         = er;
	f_cross_grad(:,pairs_range,:) = fxg;	
	
	if npairs > 20000,
		if ~rem(last_pair, chunk_size_to_use*6),
			fprintf(1, '%.0f%%', last_pair*100/npairs);
		else
			fprintf(1, '.');
		end
	end
end

clear i_x i_y i_z i_obs_x i_obs_y i_obs_z;

if npairs > 20000,
	fprintf(1, '\n');
end

% The output integrals for faces/face centers.
[m, n] = ndgrid(uint32(1:ntri), uint32(1:ntri));
tri_e_over_r = e_over_r(sub2ind([ntri ntri],m,n));
clear m n;

% Here we start populating the output arrays associated with edges. The basis
% of all arays is NxN matrix where N is the number of edges. Each entry of the
% arrays corresponds to a pair of edges: edge m which is testing
% (or observation) edge and edge n which is source edge.

% Index matrices for rows and columns. midx is the row number and therefore
% the testing edge number, nidx is the column number and therefore the source
% edge number.
[midx, nidx] = ndgrid(uint32(1:nedges), uint32(1:nedges));

% Each edge has two faces corresponding to it - positive and negative ones.
% The arrays Tm and Tn give the triangle indices for the intergal arrays
% element with indices m and n. Both are nedges-by-nedges-by-2 arrays, first two
% indices are the edge indices, third one is for positive/negative faces.
sub_m = repmat(midx, [ 1 1 2 ]);
sub_pos_neg = repmat(uint32(cat(3, 1, 2)), nedges, nedges);
Tm = uint32(edge_tris(sub2ind(size(edge_tris), sub_m, sub_pos_neg)));
sub_n = repmat(nidx, [ 1 1 2 ]);
Tn = uint32(edge_tris(sub2ind(size(edge_tris), sub_n, sub_pos_neg)));
clear midx sub_m; 

% Subscripts to build linear index of the integrals for the pairs of edges.
sub_Tm = repmat(Tm, [ 1 1 1 2 ]);
sub_Tn = repmat(permute(Tn, [ 1 2 4 3 ]), [ 1 1 2 1 ]);
clear Tm Tn; 

% Indices of the pairs of faces in the integrals evaluator output arrays. 
% nedges-by-nedges-by-2-by-2 array, first two indices are the edge indices,
% third one is the m edge triangle sign, fourth one is the m edge triangle
% sign.
iidx = sub2ind([ntri ntri], sub_Tm, sub_Tn);
clear sub_Tm;

% The e_over_r given by the integrals evaluator needs to be multiplied by
% div(f) value because the evaluator omits it.
div_f_sign = repmat(cat(4, 1, -1), [ nedges nedges 2 ]);
div_f = div_f_sign.*edge_l(repmat(nidx, [ 1 1 2 2 ]))./tri_a(sub_Tn);

% Now we can compute e_over_r for the edge pairs. Sum the integrals over the
% two faces of the source edge.
edge_e_over_r = sum(e_over_r(iidx).*div_f, 4);
clear e_over_r div_f_sign div_f;

% f_e_over_r and f_cross_grad have additional index for x-y-z components
% (both are vectors). We can not use iidx as a linear index to address these
% arrays, but we can use it when building the new linear index. Modify it
% accordingly, now it is nedges-by-nedges-by-2-by-2-by-3 array, the new
% dimension is for x-y-z.
iidx = repmat(iidx,[ 1 1 1 1 3 ]);

% x-y-z components index
xyz_idx = repmat(uint32(cat(5, 1, 2, 3)), [ nedges nedges 2 2 ]);

% f_e_over_r and f_cross_grad have additional index for the free vertex;
% vidx is a local (within a triangle) indices of the free vertex for the
% source edge n corresponding to the arrays element (m,n).
% nedges-by-nedges-by-2 array, first two indices are the edge indices,
% third one is 1 for positive and 2 for negative faces.
vidx = free_vert_loc(sub2ind(size(free_vert_loc), sub_n, sub_pos_neg));
clear sub_n sub_pos_neg;

% Permute and repmat vidx so it matches iidx and xyz_idx. The free vertex
% is one of the source triangle.
vidx = repmat(permute(uint32(vidx), [ 1 2 4 3 ]), [ 1 1 2 1 3 ]);

% Now we can compose linear index for f_e_over_r and f_cross_grad arrays.
i2idx = sub2ind(size(f_e_over_r), vidx, iidx, xyz_idx);
clear vidx iidx xyz_idx;

% f_e_over_r and f_cross_grad need to be multiplied by a normalization
% factor which is ln/(2*An) because this factor is omitted by the integrals
% evaluator. In addition, the integral over the negative face needs to be
% multiplied by -1 because the integrals evaluator directs f from free vertex
% to the evaluation point - but actually in the negative triangle it is
% pointing backwards.
f_norm_sign = repmat(cat(4, 1, -1), [ nedges nedges 2 1 3 ]);
f_norm_edge_l = edge_l(repmat(nidx, [ 1 1 2 2 3 ]));
f_norm_tri_a = tri_a(repmat(sub_Tn, [ 1 1 1 1 3 ]));
f_norm = f_norm_sign.*f_norm_edge_l./(2*f_norm_tri_a);
clear f_norm_sign f_norm_edge_l f_norm_tri_a nidx sub_Tn;

% Finally, sum the integrals over the positive and negative faces of the
% source edge after multiplying by a normalization factor
edge_f_e_over_r = sum(f_e_over_r(i2idx).*f_norm, 4);
edge_f_e_over_r = permute(edge_f_e_over_r, [ 1 2 3 5 4 ]);

% All the same like for the edge_e_over_r. We need the derivative with
% respect to the integration coordinate r, integrals evaluator treats R
% as the distance from observation to the integration points (R=|r-r'|,
% r - integration, r' - observation) so we do not need to apply minus
% sign here.
edge_f_cross_grad = sum(f_cross_grad(i2idx).*f_norm, 4);
edge_f_cross_grad = permute(edge_f_cross_grad, [ 1 2 3 5 4 ]);
