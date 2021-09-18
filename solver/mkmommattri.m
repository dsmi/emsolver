function v = mkmommattri(mesh, integf, mqo, mt, nt, integpf)
% v = mkmommattri(mesh, integf, mqo, mt, nt, integpf)
%
% Computes the moment matrix. Each element of the matrix m(m,n) is the
% corresponding operator (as defined by integf) applied to the n-th expansion
% basis function tested with respect to the m-th testing basis function.
% The expansion and testing functions are piecewise constant ones associated
% with triangles. Notice that the result omits 1/Am factor.
% The optional integpf can be used to post-process the results of integf, it
% receives the result of integf and the source and observation triangle
% indices.
%

% Number of source and observation triangles.
nmt = length(mt);
nnt = length(nt);

% Default post-integration function, does nothing
if ~exist('integpf')
	integpf = @( v, srct, obst ) ( v );
end

% If the matrix is big enough, it can not be computed at once because
% of the memory limiations - it is computed by blocks instead.
maxbls=200; % Maximum size of the block
if nmt*nnt>maxbls*maxbls,

	fprintf(1, 'Populating moment matrix');	
	
	% Pre-allocate the overall resulting matrix
	v=zeros(nmt,nnt);
	mend = 0;
	while mend < nmt,
		mstart = mend + 1;
		mend = min(mend + maxbls, nmt);
		bmt = mt(mstart:mend);

		nend = 0;
		while nend < nnt,
			nstart = nend + 1;
			nend = min(nend + maxbls, nnt);
			bnt = nt(nstart:nend);
			
			% Calculate the subblock
			vb = mkmommattri(mesh, integf, mqo, bmt, bnt, integpf);
			
			% And put it into the resulting matrix
			v(mstart:mend,nstart:nend)=vb;
			
			%fprintf(1, '.');
		end
		
		fprintf(1, '.%.0f%%', mend*100/nmt);
	end

	fprintf(1, '\n');

	return;
end

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
idx = mesh.tri(mt,:);
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
idx = mesh.tri(nt,:);
r = cat(3, mesh.x(idx), mesh.y(idx), mesh.z(idx));

% Now it is of size [ nmt nnt nq 3 3 ]
r = repmat(permute(r, [ 4 1 5 2 3 ]), [ nmt 1 nq ]);

% Inner integral of the product
% array of size [ nmt nnt nq ]
rr = reshape(r, [], 3, 3);
robsr = reshape(robs, [], 3);
obst = reshape( repmat( reshape( mt, [], 1  ), 1,   nnt, nq ), [], 1 );
srct = reshape( repmat( reshape( nt,  1, [] ), nmt, 1  , nq ), [], 1 );
intgr = integpf( integf( rr, robsr ), srct, obst );
intg = reshape(intgr, nmt, nnt, nq);

% Repmat and reshape qw so it has the same dimensions as intg
tqw = repmat(permute(qw, [ 2 3 1 ]), [ nmt nnt 1 ]);

% Evaluate the quadrature - gives the final result. The factor 1/A is omitted.
v = sum(intg.*tqw, 3)*2;
