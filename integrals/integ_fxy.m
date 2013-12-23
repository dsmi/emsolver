function fxy = integ_fxy(k, r, robs, qN, aux)
% function fxy = integ_fxy(k, r, robs)
% 
% Calculates integrals of the div(fxy*exp(-j*k*R)/R), where
% fxy = f-z(f*z) (f without z component), f is the vector between free
% vertex and integration point: f = r-r(l) and R is the distance between
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
%   fxy  - the resulting integrals. Two-dimensional array of size N-by-3,
%          first index is the source-observation pair index, second index
%          is the local index of the free vertex.
%

% To evaluate the integral, it is divided into two parts. The divergence can
% be expanded into a sum of the surface divergence and the normal derivative
% dotted with the normal vector:
%  div(f) = div_s(f)+dot(n,d(f)/dn)
% The surface divergence part of the target integral is evaluated via
% the divergence theorem, the 'normal' part is evaluated using the
% dedicated function integ_fdn.
%

% Number of the quadrature points to use
if nargin < 4
	qN = 3;
end

% Auxiliaries, if not supplied by caller
if nargin < 5
	aux = integ_aux(r, robs, qN);
end

integrand4 = exp(-j*k*aux.qR)./aux.qR;
integral4 = sum(integrand4.*aux.qW,3).*(aux.s./2);
fxy_r123 = sum(aux.nu.*repmat(integral4, [ 1 1 3 ]), 2);

% This is the part of the integral which depends on the free vertex
[ lx ly lz ] = uncat(3, r.*repmat(fxy_r123, [ 1 3 ]));

[ qrx qry qrz ] = uncat(4, permute(aux.qr, [ 1 2 4 3 ]));

integrand4x = qrx.*integrand4;
integrand4y = qry.*integrand4;
%integrand4z = qrz.*integrand4;
integral4x = sum(integrand4x.*aux.qW,3).*(aux.s./2);
integral4y = sum(integrand4y.*aux.qW,3).*(aux.s./2);
%integral4z = sum(integrand4z.*aux.qW,3).*(aux.s./2);

[ nux nuy nuz ] = uncat(3, aux.nu);

% This is the part of the integral which depends on the observation point only
rx = sum(nux.*integral4x, 2);
ry = sum(nuy.*integral4y, 2);
%rz = sum(nuz.*integral4z, 2);

fxy = (rx(:,ones(1,3))-lx) + (ry(:,ones(1,3))-ly);% + (rz(:,ones(1,3))-lz);

% Calculate 'normal' term of the divergence
fdn = integ_fdn(k, r, robs, qN, aux);
fdn(:,:,3) = 0*fdn(:,:,3); % Drop z component
nddn = sum(fdn.*repmat(aux.n, 1, 3), 3);

% Add together
fxy = fxy+nddn;
