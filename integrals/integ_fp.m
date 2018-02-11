function fp = integ_fp(k, r, robs, qN, aux)
% function fxg = integ_fp(k, r, robs, qN, aux)
% 
% Calculates integrals of the f*exp(-j*k*R)/R, where f is the vector
% between free vertex and integration point: f = r-r(l), and R is the
% distance between observation point and integration.
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
%   fp   - the resulting integrals. Three-dimensional array of size
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
% f*p = f*exp(-j*k*R)/R), where f is the RWG basis function given by
% f = r - r(l), and r is the evaluation point coordinates and r(l)
% is the free vertex.
%
% To evaluate the integral, we divide it into two parts:
%   x1 = r*exp(-j*k*R)/R
%   x2 = r(l)*exp(-j*k*R)/R
% The first part includes the observation coordinate only, while the second
% depends on the free vertex position as well.
%
% First part - depends on the observation coordinate only. In the coordinate
% system with the X and Y axes lying in the plane of triangle it can be
% written as:
%   fp_x = x*exp(-j*k*R)
%   fp_y = y*exp(-j*k*R)
% Derivatives of -1/(j*k)*exp(-j*k*R) with respect to x and y can be seen
% to be equal to:
%   d(-1/(j*k)*exp(-j*k*R))/dx = x*exp(-j*k*R)/R
%   d(-1/(j*k)*exp(-j*k*R))/dy = y*exp(-j*k*R)/R
% This is the target funtion, so we need to integrate the
% -1/(j*k)*exp(-j*k*R) along the triangle edges to obtain the target
% surface integral value.

% To handle k=0 case
if k == 0
    exp_jk = @(r) r;
    invjk = 1;
else
    exp_jk = @(r) exp(-j*k*r);
    invjk = -1/(j*k);
end

integrand1 = exp_jk(aux.qR);
integral1 = sum(integrand1.*aux.qW,3).*(aux.s./2);
fp_obs = invjk*sum(aux.nu.*repmat(integral1, [ 1 1 3 ]), 2);

% Second part - includes the free vertex position. Can be written as:
%   fp_x = x(l)*exp(-j*k*R)/R
%   fp_y = y(l)*exp(-j*k*R)/R
% where x(l) and y(l) are coordinates of the free vertex. These are
% constants and can be taken out of the integral. Thus, we only need to
% calculate integral of:
%   f = exp(-j*k*R)/R
% and multiply it by r(l) vector.
intgp = integ_p(k, r, robs, qN, aux);
fp_r123 = (r - aux.rho_r).*repmat(intgp, [ 1 3 3 ]);

% We are done! - just sum the two parts.
fp = repmat(fp_obs, 1, 3) - fp_r123;
