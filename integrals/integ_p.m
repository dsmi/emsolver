function p = integ_p(k, r, robs, qN, aux)
% function p = integ_p(k, r, robs, qN, aux)
% 
% Calculates integrals of the exp(-j*k*R)/R, where R is the distance
% between observation point and integration.
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
%   p    - the resulting integrals. Column vector.
%

% Number of the quadrature points to use
if nargin < 4
	qN = 3;
end

% Auxiliaries, if not supplied by caller
if nargin < 5
	aux = integ_aux(r, robs, qN);
end

% The integral we are going to evaluate is:
%   f = exp(-j*k*R)/R
% Consider the vector function with the following components:
%   f_x = -1/(j*k)*x*exp(-j*k*R)/(x^2+y^2)
%   f_y = -1/(j*k)*y*exp(-j*k*R)/(x^2+y^2)
% Its divergence is equal to the target function given above, thus,
% according to the divergence theorem we can integrate this vector
% function over the triangle edges instead of evaluating the surface
% integral.
% The vector function with terms f_x and f_y given above is singular
% when rho=0; to address that we add the following correction terms:
%   g_x = 1/(j*k)*x*exp(-j*k*abs(h))/(x^2+y^2)
%   g_y = 1/(j*k)*y*exp(-j*k*abs(h))/(x^2+y^2)
% so:
%   f_x_nonsing = f_x + g_x
%   f_y_nonsing = f_y + g_y
% Now it can be seen that both f_x_nonsing and f_y_nonsing are non-singular
% i.e. they converge to a finite limit when rho -> 0.
% This does not affect the integration results because divergence of vector
% function (g_x, g_y) is zero.
% Next, we simplify the things even more. For each of the edges, we need to
% integrate the dot product of the function we just found and the edge
% normal. For each of the edges we can translate the coordinates such that
% Y axis is parallel to the edge and X is perpendicular. Thus, we only need
% to integrate X component, because in this new coordinate system the edge
% normal has no Y component. The X component is constant and is equal to
% the distance from the observation point to the line the edge lies on
% (r2_minus_rho_dot_nu).

% To handle k=0 case
if k == 0
    exp_jk = @(r) r;
    invjk = 1;
else
    exp_jk = @(r) exp(-j*k*r);
    invjk = -1/(j*k);
end

inv_q_rho_2 = 1./aux.q_rho_2;
integranda = exp_jk(aux.qR).*inv_q_rho_2;
integrandb = -inv_q_rho_2;
integrala = sum(aux.rho2line.*sum(integranda.*aux.qW,3).*(aux.s/2), 2);
m = aux.cin.*exp_jk(aux.abs_h);
integralb = m.*sum(aux.rho2line.*sum(integrandb.*aux.qW,3).*(aux.s/2), 2);
p = invjk*(integrala+integralb);
