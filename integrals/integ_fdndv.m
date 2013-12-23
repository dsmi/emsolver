function fdndv = integ_fdndv(k, r, robs, v, qN, aux)
% function fdndv = integ_fdndv(k, r, robs, v, qN, aux)
% 
% Calculates integrals of the f*d2(exp(-j*k*R)/R)/dn/dv, where f is the
% vector between free vertex and integration point: f = r-r(l), R is
% the distance between observation and integration points, n is the
% triangle normal and v is the user-supplied derivative direction.
%
% Inputs:
%   k    - wavenumber, scalar
%   r    - vertices of the triangles, array of size N-by-3-by-3, first index
%          is the index of the source-observation pair, second index is the
%          triangle vertex index, third is X-Y-Z.
%   robs - observation points, array of size N-by-3, first index is 
%          the index of the source-observation pair, second index is X-Y-Z.
%   v    - a vector of three elements, the derivative direction.
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
if nargin < 5
	qN = 3;
end

% Auxiliaries, if not supplied by caller
if nargin < 6
	aux = integ_aux(r, robs, qN);
end

%
% Very brief expalation: all the expressions below can be obtained by
% differentiation of the expressions in integ_fdn with respect to v.
%

% num_of_tri-by-1-by-3 array
vr = repmat(reshape(v, [ 1 1 3 ]), size(r, 1), 1);

% Here is the part which depends on the free vertex position.

% integrals of p/|rho|^2 and of d|h|/dn*rho/|rho|^2
p = exp(-j*k*aux.qR)./aux.qR;
inv_q_rho_2 = 1./aux.q_rho_2;
integrand1a = p.*inv_q_rho_2;
integrand1b = -inv_q_rho_2;
integral1a = sum(integrand1a.*aux.qW,3).*(aux.s/2);
%integral1a = aux.h.*sum(aux.rho2line.*integral1a, 2);
exph = exp(-j*k*aux.abs_h);
integral1b = sum(integrand1b.*aux.qW,3).*(aux.s/2);
%integral1b = exph.*aux.d_abs_h_dn.*sum(aux.rho2line.*integral1b, 2);

% derivative with respect to r
dpdr = -p.*(1+j*k*aux.qR)./aux.qR;

qrr = aux.qr - repmat(permute(robs, [ 1 3 2 ]),[ 1 3 1 qN ]);
rv = permute(sum(qrr.*repmat(vr, [ 1 3 1 qN ]), 3), [ 1 2 4 3 ]);
drdv = rv./aux.qR;
dhdv = sum(aux.n.*vr, 3);

drhodv = vr-aux.n.*repmat(dhdv, [ 1 1 3 ]);
drho2dv = permute(2*sum(aux.q_rho.*repmat(drhodv, [ 1 3 1 qN ]), 3), [ 1 2 4 3 ]);

dpdv = dpdr.*drdv;

% integral of dh/dv*rho*p/|rho|^2
integral2 = integral1a;
integral2 = aux.rho2line.*integral2;
integral2 = dhdv.*sum(integral2, 2);

% integral of h*drho/dv*p/|rho|^2
drhodv_dot_nu = sum(repmat(drhodv, 1, 3).*aux.nu, 3);
integral3 = aux.h_r.*drhodv_dot_nu.*integral1a;
integral3 = sum(integral3, 2);

% d(p/|rho|^2)/dv
integrand4 = -drho2dv.*p.*inv_q_rho_2./aux.q_rho_2+dpdv./aux.q_rho_2;
integral4 = sum(integrand4.*aux.qW,3).*(aux.s/2);
integral4 = aux.h_r.*integral4;
integral4 = aux.rho2line.*integral4;
integral4 = sum(integral4, 2);

integrand5 = inv_q_rho_2;
integral5 = drhodv_dot_nu.*sum(integrand5.*aux.qW,3).*(aux.s/2);
integral5 = sum(integral5, 2);

integrand6a = -drho2dv.*inv_q_rho_2.*inv_q_rho_2;
integrand6b = inv_q_rho_2;
integral6a = aux.rho2line.*sum(integrand6a.*aux.qW,3).*(aux.s/2);
integral6b = aux.rho2line.*sum(integrand6b.*aux.qW,3).*(aux.s/2);
integral6 = sum(integral6a, 2) - j*k*aux.d_abs_h_dn.*dhdv.*sum(integral6b, 2);

fdndv_r123 = integral2+integral3+integral4-exph.*aux.d_abs_h_dn.*(integral5+integral6);
fdndv_r123 = (r - aux.rho_r).*repmat(fdndv_r123, [ 1 3 3 ]);

integral7a = aux.h.*sum(aux.rho2line.*integral1a, 2);
integral7b = exph.*aux.d_abs_h_dn.*sum(aux.rho2line.*integral1b, 2);
integral7 = integral7a + integral7b;
fdndv_r123 = fdndv_r123-repmat(-drhodv,1,3).*repmat(integral7, [ 1 3 3 ]);

% The part which depends on the observation coordinate only
integrand8 = repmat(aux.h_r, [ 1 1 qN ]).*dpdv+p.*repmat(dhdv, [ 1 3 qN ]);
integral8 = sum(integrand8.*aux.qW,3).*(aux.s./2);
fdndv_obs = sum(aux.nu.*repmat(integral8, [ 1 1 3 ]), 2);

fdndv = repmat(fdndv_obs, 1, 3) - fdndv_r123;
