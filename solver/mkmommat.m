function v = mkmommat(mesh, integf, mqo, me, ne, fprod)
% v = mkmommat(mesh, integf, mqo, me, ne, fprod)
%
% Computes the moment matrix. Each element of the matrix m(m,n) is the
% corresponding operator (as defined by integf) applied to the n-th expansion
% basis function tested with respect to the m-th testing basis function.
% The expansion and testing functions are the same and are RWG, which are
% the linear basis functions associated with edges and defined over two
% triangles sharing a common edge.
% This function allows to pass handle of a funcion which evaluates
% the testing product. It is given the tested value and testing function
% which are arrays of size [ M 2 N 2 nq 3 ] and is expected to return
% the array of size [ M 2 N 2 nq ] with the testing product value.
% By default the dot product is calculated.
%

% Function which evaluates the testing product, the default value
if ~exist('fprod')
	fprod = @(f,g) sum(f.*g, 6);
end

M = length(me);
N = length(ne);

% If the matrix is big enough, it can not be computed at once because
% of the memory limiations - it is computed by blocks instead.
maxbls=200; % Maximum size of the block
if M*N>maxbls*maxbls,

	fprintf(1, 'Populating moment matrix');	
	
	% Pre-allocate the overall resulting matrix
	v=zeros(M,N);
	mend = 0;
	while mend < M,
		mstart = mend + 1;
		mend = min(mend + maxbls, M);
		bme = me(mstart:mend);

		nend = 0;
		while nend < N,
			nstart = nend + 1;
			nend = min(nend + maxbls, N);
			bne = ne(nstart:nend);
			
			% Calculate the subblock
			vb = mkmommat(mesh, integf, mqo, bme, bne);
			
			% And put it into the resulting matrix
			v(mstart:mend,nstart:nend)=vb;
			
			%fprintf(1, '.');
		end
		
		fprintf(1, '.%.0f%%', mend*100/M);
	end

	fprintf(1, '\n');

	return;
end

% Triangles referred by m (testing) and n (source) edges. M-by-2 and N-by-2
% arrays correspondingly. 
mt = mesh.edge_tris(me,:);
nt = mesh.edge_tris(ne,:);

% Column vectors mtu and ntu define local subsets of triangles referred
% by m and n edges.
mtu = unique(mt(:));
ntu = unique(nt(:));

% Number of triangles referred by m and n edges
nmt = length(mtu);
nnt = length(ntu);

% Given the triangle index in the full triangles vector mt_a2l and
% nt_a2l vectors allow to find index of the corresponding triangle in
% mtu and ntu vector (subsets of triangles referred by m and n edges).
% If the triangle is not referred, the index is zero.
mt_a2l = zeros(size(mesh.tri,1),1);
mt_a2l(mtu) = 1:length(mtu);
nt_a2l = zeros(size(mesh.tri,1),1);
nt_a2l(ntu) = 1:length(ntu);

% Triangles referred by m (testing) and n (source) edges, but this time
% the index refers to the local subset of the triangles, see mtu and ntu.
mtl = reshape(mt_a2l(mt), size(mt)); % the reshape is needed if M=1
ntl = reshape(nt_a2l(nt), size(nt));

% Quadrature points
qo = mqo; % Order of the quadrature
if qo < 2,
	nq = 1;
	qa = [ 1/3 1/3 ];
	qw = 0.5;
else
	nq = qo.*qo;
	[qa,qw]=simplexquad(qo,2);
endif

% This gives us barycentirc coordinates, now qa is nq-by-3 array.
qa=[ qa 1-sum(qa,2) ];

% Observation points for the integrals evaluator.
qa = repmat(permute(qa,[ 3 2 1 ]), nmt, 1);
idx = mesh.tri(mtu,:);
qx = sum(repmat(mesh.x(idx), [ 1 1 nq ]).*qa, 2);
qy = sum(repmat(mesh.y(idx), [ 1 1 nq ]).*qa, 2);
qz = sum(repmat(mesh.z(idx), [ 1 1 nq ]).*qa, 2);

% Quadrature points, size is [ nmt 3 nq ], second index is X-Y-Z.
qr = [ qx qy qz ];

% Observation points for the integrals evaluator, size is [ nmt nnt nq 3 ].
robs = repmat(permute(qr, [ 1 4 3 2 ]), 1, nnt);

% Source triangles for the integrals evaluator.
% array of size [ nnt 3 3 ], first index is the triangle index,
% second is the vertex, third is X-Y-Z.
idx = mesh.tri(ntu,:);
r = cat(3, mesh.x(idx), mesh.y(idx), mesh.z(idx));

% Now it is of size [ nmt nnt nq 3 3 ]
r = repmat(permute(r, [ 4 1 5 2 3 ]), [ nmt 1 nq ]);

% Inner integral of the product
% array of size [ nmt nnt nq 3 3 ]. The fourth index is the free vertex
% index, the fifth is X-Y-Z.
rr = reshape(r, [], 3, 3);
robsr = reshape(robs, [], 3);
intgr = integf(rr, robsr);
intg = reshape(intgr, nmt, nnt, nq, 3, 3);

% Now build array of integrals associated with edge triangles. Its size is
% [ M 2 N 2 nq 3 ]. First index is the testing edge index,
% second index is the source edge index, third index the pos/neg triangle
% of the testing edge, fourth one is pos/neg triangle of the source edge,
% fifth is the quadrature point index, sixth is X-Y-Z.
submt = repmat(mtl, [ 1 1 N 2 nq 3 ]);
subnt = repmat(permute(ntl, [ 3 4 1 2 ]), [ M 2 1 1 nq 3 ]);
subq = repmat(permute(1:nq, [ 1 3 4 5 2 ]), [ M 2 N 2 1 3 ]);
subf = repmat(permute(mesh.free_vert_loc(ne, :), [ 3 4 1 2 ]), [ M 2 1 1 nq 3 ]);
subx = repmat(permute(1:3, [ 1 3 4 5 6 2 ]), [ M 2 N 2 nq ]);
intgn = intg(sub2ind(size(intg), submt, subnt, subq, subf, subx));

% Quadrature points associated with positive and negative triangles
% of edge m. size is [ M 2 nq 3 ]; first index is edge, second is pos/neg
% triangle, third is the quadrature index, fourth is xyz.
subt = repmat(mtl, [ 1 1 nq 3 ]);
subq = repmat(permute(1:nq, [ 1 3 2 4 ]), [ M 2 1 3 ]);
subx = repmat(permute(1:3, [ 1 3 4 2 ]), [ M 2 nq ]);
qrm = qr(sub2ind(size(qr), subt, subx, subq));

% Free vertices associated with positive and negative triangles
% of edge m. size is [ M 2 3 ]; first index is edge, second is pos/neg
% triangle, third is xyz.
idx = mesh.free_vert(me,:);
fvrm = cat(3, mesh.x(idx), mesh.y(idx), mesh.z(idx));

% Permute and repmat it so it has the same dimensions as qrm
fvrm = repmat(permute(fvrm, [ 1 2 4 3 ]), [ 1 1 nq ]);

% Vector from free vertex to the quadrature points (and in negative direction
% for the negative triangle). [ M 2 nq 3 ]; first index is edge, second
% is pos/neg triangle, third is the quadrature index, fourth is xyz.
rhom = (qrm-fvrm).*repmat([ 1 -1 ], [ M 1 nq 3]);

% Permute and repmat it so it has the same dimensions as intgn
rhom = repmat(permute(rhom, [ 1 2 5 6 3 4 ]), [ 1 1 N 2 ]);

% Evaluate the testing product
tprodq = fprod(rhom, intgn);

% Repmat and reshape qw so it has the same dimensions as tprodq
tqw = repmat(permute(qw, [ 2 3 4 5 1 ]), [ M 2 N 2 1 ]);

% Evaluate the quadrature
tprod = sum(tprodq.*tqw, 5);

% The integrals needs to be multiplied by normalization factor which is
% ln/(2*An(+/-)) because this factor is omitted by the integrals evaluator.
% In addition, the integral over the negative face needs to be multiplied
% by -1 because the integrals evaluator directs f from free vertex to
% the evaluation point - but actually in the negative triangle it is
% pointing backwards.
fn_sign = repmat(cat(4, 1, -1), [ M 2 N ]);
fn_edge_l = repmat(permute(mesh.edge_l(ne), [ 2 3 1 ]), [ M 2 1 2 ]);
fn_tri_a = repmat(permute(mesh.edge_tri_a(ne, :), [ 3 4 1 2 ]), M, 2);
fn_norm = fn_sign.*fn_edge_l./(2*fn_tri_a);

% Normalize and sum values for pos and neg triangles of edge n
tprodn = sum(tprod.*fn_norm, 4);

% One more normalization, now for the testing edge. There is no need to
% multiply by 1/(2*Am) because it was omitted when evaluating the
% quadrature.
fm_edge_l = repmat(mesh.edge_l(me), [ 1 2 N ]);
fm_norm = fm_edge_l;

% Sum values for pos and neg triangles of edge m
tprodmn = sum(tprodn.*fm_norm, 2);

v = permute(tprodmn, [ 1 3 2 ]);
