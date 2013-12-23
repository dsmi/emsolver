function fdn = integ_fdn(k, r, robs, qN, aux)
% function fdn = integ_fdn(k, r, robs)
% 
% Calculates integrals of the f*d(exp(-j*k*R)/R)/dn, where f is the
% vector between free vertex and integration point: f = r-r(l) and R is
% the distance between observation and integration points.
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
%   fdn  - the resulting integrals. Three-dimensional array of size
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
%  fdn = f*d(exp(-j*k*R)/R)/dn, where f is the RWG basis function given
%    by f = r - r(l), r is the evaluation point coordinates, r(l) is
%    the free vertex and n is the normal to the plane of triangle
% The integral is divided into two parts:
%   x1 = r*d(exp(-j*k*R)/R)/dn
%   x2 = r(l)*d(exp(-j*k*R)/R)/dn
% The first part includes the observation coordinate only, while the second
% depends on the free vertex position as well.
% To address the first part:
%  div_s(h*exp(-j*k*R)/R) = rho*d(exp(-j*k*R)/R)/dn, and the surface integral
%    is transformed to the line integral via the divergence theorem.
% The second part, which depends on the free vertex position:
%  div_s((rho/rho^2)*h*exp(-j*k*R)/R) = d(exp(-j*k*R)/R)/dn
% Similar, the divergence theorem is used. To get rid of the singularity
% the following correction term is added to the integrand
%   g = (rho/rho^2)*h*exp(-j*k*|h|)/|h|, where h is the distance from the
% observation point to the plane of triangle
% 

% Here is the second part, which depends on the free vertex position.
p = exp(-j*k*aux.qR)./aux.qR;
inv_q_rho_2 = 1./aux.q_rho_2;
integrand4a = p.*inv_q_rho_2;
integrand4b = -inv_q_rho_2;
integral4a = aux.h.*sum(aux.rho2line.*sum(integrand4a.*aux.qW,3).*(aux.s/2), 2);
m = aux.d_abs_h_dn.*exp(-j*k*aux.abs_h);
integral4b = m.*sum(aux.rho2line.*sum(integrand4b.*aux.qW,3).*(aux.s/2), 2);
fdn_r123 = (r - aux.rho_r).*repmat(integral4a+integral4b, [ 1 3 3 ]);

% The part which depends on the observation coordinate only
integrand5 = p;
integral5 = aux.h_r.*sum(integrand5.*aux.qW,3).*(aux.s./2);
fdn_obs = sum(aux.nu.*repmat(integral5, [ 1 1 3 ]), 2);

fdn = repmat(fdn_obs, 1, 3) - fdn_r123;
