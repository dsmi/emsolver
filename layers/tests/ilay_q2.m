function [ fxg ] = ilay_q2(lay, gi, r, robs, qN)
%
% Computes all the same intagrals as computed by ilay_fp and ilay_p
% using a very fine quadrature. Works very slow, can not handle self-terms
% properly, used to test these functions.
%

% Number of quadrature points for the line integrals over edges.
if nargin < 5
	qN = 32;
end

num_of_tri = size(r,1);

% Compute the triangle normals and areas
edge12 = r(:,2,:) - r(:,1,:);
edge23 = r(:,3,:) - r(:,2,:);
normals = cross(edge12, edge23, 3);
nl = sqrt(sum(normals.*normals,3)); % Length of the normals.
normals = normals./repmat(nl, [ 1 1 3 ]); % Normalize
tri_a = nl/2;

fxg = zeros(num_of_tri,3,3);

for itri=1:num_of_tri,

	% Vertices of the triangle.
	x = r(itri,:,1);
	y = r(itri,:,2);
	z = r(itri,:,3);

	% Observation point
	xobs = robs(itri,1);
	yobs = robs(itri,2);
	zobs = robs(itri,3);

	% Quadrature points
	[a,w]=simplexquad(qN,2);
	a = [ a 1-sum(a,2) ]; % barycentric coordinates
	qn = size(a,1);
	qX = sum(x(ones(qn,1),:).*a, 2);
	qY = sum(y(ones(qn,1),:).*a, 2);
	qZ = sum(z(ones(qn,1),:).*a, 2);
	w = w*tri_a(itri)*2;

	% Observation and source poits
	ro = [ qX.'; qY.'; qZ.' ];
	rs = [ xobs(1,ones(1,qn)); yobs(1,ones(1,qn)); zobs(1,ones(1,qn)) ];

	% Calculate integral of the scalar potential (as computed by ilay_p)
	p_qpt = dfield_v(lay,gi,ro,rs);
	p(itri)=sum(w.*shiftdim(p_qpt, 3));

	[ dFvv dFzx dFzy dFzz ] = dfield_df(lay,gi,ro,rs);
	dFvv = permute(dFvv, [ 4 1 2 3 ]);
	dFzx = permute(dFzx, [ 2 1 ]);
	dFzy = permute(dFzy, [ 2 1 ]);
	dFzz = permute(dFzz, [ 4 1 2 3 ]);

	for l=1:3,
		% Evaluate f(r) at the quadrature points
		fx=qX-x(l);
		fy=qY-y(l);
		fz=qZ-z(l);
		fxr=repmat(qX-x(l), 1, 3);
		fyr=repmat(qY-y(l), 1, 3);
		fzr=repmat(qZ-z(l), 1, 3);
		
		dffx = dFvv.*fxr;
		dffy = dFvv.*fyr;
		dffz = dFzx.*fxr+dFzy.*fyr+dFzz.*fzr;

		% curl(f*F). Minus, because curl(f*g) = -cross(f, grad(g))
		fxg(itri,l,1) = -sum((dffz(:,2)-dffy(:,3)).*w);
		fxg(itri,l,2) = -sum((dffx(:,3)-dffz(:,1)).*w);
		fxg(itri,l,3) = -sum((dffy(:,1)-dffx(:,2)).*w);
	end
end
