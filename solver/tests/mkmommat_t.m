function v = mkmommat_t(mesh, integf, mqo, me, ne)
% v = mkmommat_t(mesh, integf, mqo, me, ne)
%
% Does all the same what mkmommat does, but in a straightforward
% non-vectorized manner. Used to test mkmommat
%

M = length(me);
N = length(ne);

% Triangles referred by m (testing) and n (source) edges. M-by-2 and N-by-2
% arrays correspondingly. 
mt = mesh.edge_tris(me,:);
nt = mesh.edge_tris(ne,:);

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

tsgn = [ 1 -1 ];

result = zeros(M,N);

for im=1:M,
	for in=1:N,
		m=me(im);
		n=ne(in);
		for imt=1:2,
			for int=1:2,
				% Indices of the source and testing triangles
				mtt = mt(im,imt);
				ntt = nt(in,int);
				% Source triangle
				rx = mesh.x(mesh.tri(ntt, :));
				ry = mesh.y(mesh.tri(ntt, :));
				rz = mesh.z(mesh.tri(ntt, :));
				r = cat(3, rx, ry, rz);
				tprod = 0;
				for iq=1:nq,
					% Observation point
					ox = sum(mesh.x(mesh.tri(mtt,:)).*qa(iq,:), 2);
					oy = sum(mesh.y(mesh.tri(mtt,:)).*qa(iq,:), 2);
					oz = sum(mesh.z(mesh.tri(mtt,:)).*qa(iq,:), 2);
					robs = [ ox oy oz ];
					% Integral value
					intg = integf(r, robs);
					% Free vertex
					fx = mesh.x(mesh.free_vert(m,imt));
					fy = mesh.y(mesh.free_vert(m,imt));
					fz = mesh.z(mesh.free_vert(m,imt));
					% Vector from free vertex to observation point
					rho = cat(3, ox - fx, oy - fy, oz - fz)*tsgn(imt);
					% Testing dot product
					tpq = sum(rho.*intg(:,mesh.free_vert_loc(n,int),:),3);
					tprod = tprod + tpq*qw(iq);
				end
				fn_norm = mesh.edge_l(n)/(2*mesh.tri_a(ntt))*tsgn(int);
				fm_norm = mesh.edge_l(m);
				result(im,in) = result(im,in)+tprod*fn_norm*fm_norm;
			end
		end
	end
end

v = result;
