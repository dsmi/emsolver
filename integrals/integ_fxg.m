function fxg = integ_fxg(k, r, robs, qN, aux)
% function fxg = integ_fxg(k, r, robs)
% 
% Calculates integrals of the f X grad(exp(-j*k*R)/R), where f is the
% vector between free vertex and integration point: f = r-r(l), X means
% cross product, and R is the distance between observation and integration
% points.
%
% Inputs:
%   k    - wavenumber, scalar
%   r    - vertices of the triangles, array of size N-by-3-by-3, first index
%          is the index of the source-observation pair, second index is the
%          triangle vertex index, third is X-Y-Z.
%   robs - observation points, array of size N-by-3, first index is 
%          the index of the source-observation pair, second index is X-Y-Z.
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
if nargin < 4
	qN = 3;
end

% Auxiliaries, if not supplied by caller
if nargin < 5
	aux = integ_aux(r, robs, qN);
end

% We compute the integral of:
% f_cross_grad = f X grad(exp(-j*k*R)/R), where f is the RWG basis function
%  given by f = r - r(l), and r is the evaluation point coordinates and
%  r(l) is the free vertex.
%
% First, we compute the tangential (to the plane of triangle) part. Similar
% to the integrals of f*p, we divide it into two parts:
%   x1 = r X grad(exp(-j*k*R)/R)
%   x2 = r(l) X grad(exp(-j*k*R)/R)
% The first part includes the observation coordinate only, while the second
% depends on the free vertex position as well.
%
% First part - depends on the observation coordinate only. The tangential
% part of the target function can be written as:
%   fxg_t_x = y*d(exp(-j*k*R)/R)/dz
%   fxg_t_y = -x*d(exp(-j*k*R)/R)/dz
% Derivatives of z*exp(-j*k*R)/R by x and y can be seen to be equal to:
%   d(z*exp(-j*k*R)/R)/dx = x*d(exp(-j*k*R)/R)/dz
%   d(z*exp(-j*k*R)/R)/dy = y*d(exp(-j*k*R)/R)/dz
% Vector cross product of this with the triangle normal gives the target
% funtion, so we need to integrate the z*exp(-j*k*R)/R along the triangle
% edges to obtain the target surface integral vaule (see aforementioned
% papers on Gordon Bilow transform for additional details)
%integrand3 = exp(-j*k*aux.qR)./aux.qR;
%integral3 = sum(integrand3.*aux.qW,3).*(aux.s./2);
%fxg_obs = sum(aux.nu.*repmat(integral3, [ 1 1 3 ]), 2);
%fxg_obs = cross(, fxg_obs, 3);
% ! fxg_obs has been merged with fxg_n

% Second part - includes the free vertex position. Can be written as:
%   fxg_t_x = y(l)*d(exp(-j*k*R)/R)/dz
%   fxg_t_y = -x(l)*d(exp(-j*k*R)/R)/dz
% where x(l) and y(l) are coordinates of the free vertex. These are
% constants and can be taken out of the integral. Thus, we only need to
% calculate integral of:
%   f = d(exp(-j*k*R)/R)/dz
% over the triangle. Consider the vector function with the following
% components:
%   f_x = x*z*exp(-j*k*R)/(R*(x^2+y^2))
%   f_y = y*z*exp(-j*k*R)/(R*(x^2+y^2))
% Its divergence is equal to the target function given above, thus,
% according to the divergence theorem we can integrate this vector
% function over the triangle edges instead of evaluating the surface
% integral.
p = exp(-j*k*aux.qR)./aux.qR;
inv_q_rho_2 = 1./aux.q_rho_2;
integrand4a = p.*inv_q_rho_2;
integrand4b = -inv_q_rho_2;
integrand4 = integrand4a - integrand4b;
integral4a = aux.h.*sum(aux.rho2line.*sum(integrand4a.*aux.qW,3).*(aux.s/2), 2);
m = aux.d_abs_h_dn.*exp(-j*k*aux.abs_h);
integral4b = m.*sum(aux.rho2line.*sum(integrand4b.*aux.qW,3).*(aux.s/2), 2);
fxg_r123 = (r - aux.rho_r).*repmat(integral4a+integral4b, [ 1 3 3 ]);
fxg_r123 = cross(aux.n_r, fxg_r123, 3);

% The tangential part is done!
%f_cross_grad_t = repmat(fxg_obs, 1, 3) - fxg_r123;
% ! fxg_obs has been merged with fxg_n
f_cross_grad_t = fxg_r123;

% Now compute the normal part. Similar to all the other integrals we divide
% it into two parts, the first part includes the observation coordinate only,
% and the second depends on the free vertex position.
% The first part is given by:
%  f1 = x*dgdy-y*dgdx
% f1 is always zero, so we only need to compute the integral of the second
% part which involves the free vertex position. The rest is fairly
% obvious I believe.
integrand5 = p;
integral5 = sum(integrand5.*aux.qW,3).*(aux.s./2);
fxg_n = sum(aux.nu.*repmat(integral5, [ 1 1 3 ]), 2);
%f_cross_grad_n = -cross(r - aux.rho_r, repmat(fxg_n, 1, 3), 3);
% ! fxg_obs has been added to fxg_n

r_minus_robs = r - repmat(permute(robs, [ 1 3 2 ]),[ 1 3 1 ]);
f_cross_grad_n = -cross(r_minus_robs, repmat(fxg_n, 1, 3), 3);

fxg = f_cross_grad_t + f_cross_grad_n;
