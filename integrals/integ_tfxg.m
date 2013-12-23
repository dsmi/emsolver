function fxg = integ_tfxg(k, r, robs, M, qN, aux)
% function fxg = integ_tfxg(k, r, robs, M, qN, aux)
% 
% Calculates integrals of the 
%   M*f X grad(exp(-j*k*R)/R),
% where f is the vector between free vertex and integration point: f = r-r(l),
% X means cross product, R is the distance between observation and integration
% points, and M is the caller-provided transform matrix, which is applied to f.
% If M is the identity matrix the result of this function is the same as
% integ_fxg. The current implementation takes into accout only the diagonal
% elements of M, so the only possible transform is the scale transform.
%
% Inputs:
%   k    - wavenumber, scalar
%   r    - vertices of the triangles, array of size N-by-3-by-3, first index
%          is the index of the source-observation pair, second index is the
%          triangle vertex index, third is X-Y-Z.
%   robs - observation points, array of size N-by-3, first index is 
%          the index of the source-observation pair, second index is X-Y-Z.
%   M    - 3-by-3 transform matrix to be applied to f. The current
%          implementation takes into accout only the diagonal elements of M.
%   qN   - The target surface integral is evaluated by reduction to the
%          line integrals over edges. qN is the number of quadrature points
%          to use over each edge. This parameter is optional, the default
%          value is 3.
%   aux  - structure carrying auxiliary values, as returned by integ_aux.
%          When this function is called multiple times with the same r
%          and robs but with different k, it may be reasonable from
%          the performance point of view to compute the auxiliary values
%          once and use them for all the calls.
%          This parameter is optional, if aux is not supplied by the caller
%          integ_aux is called.
%   
% Outputs:
%   fxg  - the resulting integrals. Three-dimensional array of size
%          N-by-3-by-3, first index is the source-observation pair index,
%          second index is the local index of the free vertex,
%          third is X-Y-Z.
%

% Number of the quadrature points to use
if nargin < 5
	qN = 3;
end

% Auxiliaries, if not supplied by caller
if nargin < 6
	aux = integ_aux(r, robs, qN);
end

% We compute the integral of:
% f_cross_grad = M*f X grad(exp(-j*k*R)/R), where f is the RWG basis function
%  given by f = r - r(l), and r is the evaluation point coordinates and
%  r(l) is the free vertex.
%
% First, we compute the integrals of fx*grad(g), fy*grad(g) and fz*grad(g), 
%  where fx = x*dot(f,x), x is the basis vector, and g = exp(-j*k*R)/R.
% It can be computed as:
%  fx*grad(g) = grad(fx*g) - g*grad(fx)
% Next:
%  grad(fx*g) = grad_s(fx*g) + d/dn(fx*g), grad_s is the surface gradient.

% Surface gradient of exp(-j*k*R)/R
integrand4 = exp(-j*k*aux.qR)./aux.qR;
integral4 = sum(integrand4.*aux.qW,3).*(aux.s./2);
gsg = sum(aux.nu.*repmat(integral4, [ 1 1 3 ]), 2);

% This is the part of the integral of grad_s(f(x/y/z)*exp(-j*k*R)/R) which
% depends on the free vertex
gsfxg_r123 = repmat(r(:,:,1), [ 1 1 3 ]).*repmat(gsg, [ 1 3 ]);
gsfyg_r123 = repmat(r(:,:,2), [ 1 1 3 ]).*repmat(gsg, [ 1 3 ]);
gsfzg_r123 = repmat(r(:,:,3), [ 1 1 3 ]).*repmat(gsg, [ 1 3 ]);

[ qrx qry qrz ] = uncat(4, permute(aux.qr, [ 1 2 4 3 ]));
integrand4x = qrx.*integrand4;
integrand4y = qry.*integrand4;
integrand4z = qrz.*integrand4;
integral4x = sum(integrand4x.*aux.qW,3).*(aux.s./2);
integral4y = sum(integrand4y.*aux.qW,3).*(aux.s./2);
integral4z = sum(integrand4z.*aux.qW,3).*(aux.s./2);

% This is the part of the integral of grad_s(f(x/y/z)*exp(-j*k*R)/R)
% which depends on the observation point only
gsfxg_r = sum(aux.nu.*repmat(integral4x, [ 1 1 3 ]), 2);
gsfyg_r = sum(aux.nu.*repmat(integral4y, [ 1 1 3 ]), 2);
gsfzg_r = sum(aux.nu.*repmat(integral4z, [ 1 1 3 ]), 2);

% This is the integral of grad_s(f(x/y/z)*g)
gsfxg = gsfxg_r(:,ones(1,3),:) - gsfxg_r123;
gsfyg = gsfyg_r(:,ones(1,3),:) - gsfyg_r123;
gsfzg = gsfzg_r(:,ones(1,3),:) - gsfzg_r123;

% Integrals of d/dn(fx*g), d/dn(fy*g), d/dn(fz*g)
fdn = integ_fdn(k, r, robs, qN, aux);
ddnfxg = fdn(:,:,[ 1 1 1 ]).*aux.n_r;
ddnfyg = fdn(:,:,[ 2 2 2 ]).*aux.n_r;
ddnfzg = fdn(:,:,[ 3 3 3 ]).*aux.n_r;

% Gradients of fx, fy and fx
dfx = repmat(cat(3, 1, 0, 0), size(r, 1), 1) - aux.n(:,:,[ 1 1 1 ]).*aux.n;
dfy = repmat(cat(3, 0, 1, 0), size(r, 1), 1) - aux.n(:,:,[ 2 2 2 ]).*aux.n;
dfz = repmat(cat(3, 0, 0, 1), size(r, 1), 1) - aux.n(:,:,[ 3 3 3 ]).*aux.n;

% Integrals of p*grad(fx), p*grad(fy), p*grad(fz)
p = integ_p(k, r, robs, qN, aux);
pr = repmat(p, [ 1 1 3 ]);
gdfx = pr.*dfx;
gdfy = pr.*dfy;
gdfz = pr.*dfz;

fxdg = gsfxg - repmat(gdfx, 1, 3) + ddnfxg;
fydg = gsfyg - repmat(gdfy, 1, 3) + ddnfyg;
fzdg = gsfzg - repmat(gdfz, 1, 3) + ddnfzg;

% Finally, given the integrals of fx*grad(g), fy*grad(g) and fz*grad(g)
% compute cross(f,grad(g))
fxg = 0*fxdg;
fxg(:,:,1) = M(2,2)*fydg(:,:,3) - M(3,3)*fzdg(:,:,2);
fxg(:,:,2) = M(3,3)*fzdg(:,:,1) - M(1,1)*fxdg(:,:,3);
fxg(:,:,3) = M(1,1)*fxdg(:,:,2) - M(2,2)*fydg(:,:,1);
