function dp = integ_dp(k, r, robs, qN, aux)
% function dp = integ_dp(k, r, robs, qN, aux)
% 
% Calculates integrals of the grad(exp(-j*k*R)/R), R is the distance between
% observation and integration points.
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
%   dp   - the resulting integrals. Three-dimensional array of size
%          N-by-1-by-3, first index is the source-observation pair index,
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
% dp = grad(exp(-j*k*R)/R)
%
% First, we compute the tangential (to the plane of triangle) part, which is
% the integral of the surface gradient, via the gradient theorem.
p = exp(-j*k*aux.qR)./aux.qR;
integrand3 = p;
integral3 = sum(integrand3.*aux.qW,3).*(aux.s./2);
dpt = sum(aux.nu.*repmat(integral3, [ 1 1 3 ]), 2);

% The normal part is coputed via the divergence theorem using:
%  d(exp(-j*k*R)/R)/dn = div_s(rho/rho^2*h*exp(-j*k*R)/R)
% where the div_s is the surface divergence.
inv_q_rho_2 = 1./aux.q_rho_2;
integrand4a = p.*inv_q_rho_2;
integrand4b = -inv_q_rho_2;
integrand4 = integrand4a - integrand4b;
integral4a = aux.h.*sum(aux.rho2line.*sum(integrand4a.*aux.qW,3).*(aux.s/2), 2);
m = aux.d_abs_h_dn.*exp(-j*k*aux.abs_h);
integral4b = m.*sum(aux.rho2line.*sum(integrand4b.*aux.qW,3).*(aux.s/2), 2);
dpn = aux.n.*repmat(integral4a+integral4b, [ 1 1 3 ]);

% Sum the tangential and the normal parts
dp = dpt+dpn;
