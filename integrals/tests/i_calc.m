function [ f_e_over_r, e_over_r, f_cross_grad, one_over_r, ...
           e1_over_r, f_over_r, f_e1_over_r ] ...
               = i_calc(k, r, robs, qn)
%
% !This function is now deprecated and is used for the tests only.
% The integ_p and integ_fp are used instead.
% 
% Computes potential integrals.
%
% k : Wavenumber
%
% r : vertices of the triangles, array of size N-by-3-by-3, first index
% is the index of the source-observation pair, second index is the triangle
% vertex index, third is X-Y-Z.
% 
% robs : observation points, array of size N-by-3, first index is the index
% of the source-observation pair, second index is X-Y-Z. 
% 
% f_e_over_r : integrals of the f*exp(-j*k*R)/R, where f is the vector
% between free vertex and integration point: f = r-r(l), and R is the
% distance between observation point and integration. Three-dimensional
% array of size 3-by-num_of_triangles-by-3, first index is the local index
% of the free vertex in triangle, second is the triangle index, third
% index is X-Y-Z.
%
% e_over_r : integrals of the exp(-j*k*R)/R, where R is the distance
% between observation point and integration. Row vector of length
% num_of_triangles.
%
% f_cross_grad : integrals of the f X grad(exp(-j*k*R)/R), where f is the
% vector between free vertex and integration point: f = r-r(l), X means
% cross product, and R is the distance between observation and integration
% points. Three-dimensional array of size 3-by-num_of_triangles-by-3,
% first index is the local index of the free vertex in triangle, second 
% is the triangle index, third index is X-Y-Z.
%
% one_over_r : Auxiliary output which is used only for testing. Integrals
% of the 1/R, where R is the distance between observation point and
% integration. Row vector of length num_of_triangles.
%
% e1_over_r : auxiliary output which is used only for testing. Integrals of
% the (exp(-j*k*R)-1)/R, where R is the distance between observation
% point and integration. Row vector of length num_of_triangles.
%
% f_over_r : auxiliary output which is used only for testing. Integrals of
% f/R, where f is the vector between free vertex and integration point:
% f = r-r(l), and R is the distance between observation point and integration.
% Three-dimensional array of size 3-by-num_of_triangles-by-3, first index is
% the local index of the free vertex in triangle, second  is the triangle
% index, third index is X-Y-Z.
%
% f_e1_over_r : auxiliary data which is used only for testing. Integrals of
% f*(exp(-j*k*R)-1)/R, where f is the vector between free vertex and
% integration point: f = r-r(l), and R is the distance between observation
% point and integration. Three-dimensional array of size 
% 3-by-num_of_triangles-by-3, first index is the local index of the free
% vertex in triangle, second  is the triangle index, third index is X-Y-Z.
% 
% All the integrals are computed using Gordon-Bilow transform except
% f_cross_grad; for more details on the Gordon-Bilow transform see:
% W. B. Gordon and H. J. Bilow, "Reduction of surface integrals to
% contour integrals," IEEE Trans. Antennas Propag., vol. 50, no. 3,
% pp. 308–311, 2002.
% W. B. Gordon and H. J. Bilow, "Line-Integral Approach to Computing
% Impedance Matrix Elements" IEEE Trans. Antennas Propag., vol. 55,
% no. 10, pp. 2767–2772, October 2007.


% Number of triangles.
N = size(r, 1);

% Vertices. 3-by-N-by-3 array. First index is the local index of the vertex
% in triangle, second index is the triangle index, third index is X-Y-Z.
r = permute(r, [ 2 1 3 ]);

% Observation points. 1-by-N-by-3 array. Second index is the observation
% point index, third index is X-Y-Z.
robs = shiftdim(robs, -1);

% Edge endpoints. First edge is one opposite to the first vertex (bounded by
% vertices 2 and 3), second is one between vertices 3 and 1, third is
% between 1 and 2. 3-by-N-by-3 array. First index is the local index of the
% egde in triangle, second one is the triangle index, third index is X-Y-Z.
r1 = zeros(3, N, 3);
r1(1,:,:) = r(2,:,:);
r1(2,:,:) = r(3,:,:);
r1(3,:,:) = r(1,:,:);
r2 = zeros(3, N, 3);
r2(1,:,:) = r(3,:,:);
r2(2,:,:) = r(1,:,:);
r2(3,:,:) = r(2,:,:);

% Edge vectors. 3-by-N-by-3 array. First index is the local index of the egde
% in triangle, second one is the triangle index, third index is X-Y-Z.
edges = r2 - r1;

% Edge lengths. 3-by-N array. First index is the local index of the egde in
% triangle, second one is the triangle index.
s = sqrt(sum(edges.^2,3));

% Edge tangentials - normalized edge vectors.
t = edges ./ s(:,:,ones(1,3));

% Triangle normals, non-normalized yet. 1-by-N-by-3 array, second index is
% the triangle index, third index is X-Y-Z.
nn = cross(t(1,:,:),t(2,:,:),3);

% Length of the normals, row vector.
nl = sqrt(sum(nn.^2,3));

% Normalized normals. 1-by-N-by-3 array.
n = nn ./ nl(:,:,ones(1,3));

% Free coefficient in the plane equation of the plane of triangle. The plane
% equation is a*x + b*y + c*z + d = 0, where [ a b c ] it the normal to the
% plane and d is the free coefficient. Row vector. Is needed to compute
% projection of the observation point to the triangle plane. Row vector.
d = -sum(r(1,:,:).*n(1,:,:),3);

% Signed distance from observation point to the plane of triangle.
% Row vector.
h = sum(robs.*n,3) + d;
% h matrix repeated - 3-by-N array, N is the triangle index
h_r = h(ones(1,3),:);

% Projection of the observation point to the plane of triangle. 1-by-N-by-3
% array, second index is the triangle index, third index is X-Y-Z.
rho = robs - n .* h(:,:,ones(1,3));

% Edge normal, cross product of the edge tangent and Z vector. First
% index is the local index of the egde in triangle, second one is the
% triangle index, third index is X-Y-Z.
nu = cross(t, n(ones(1,3),:,:), 3);

% robs matrix repeated so it has the same dimensions as r, r1, r2
robs_r = robs(ones(1,3),:,:);

% rho matrix repeated so it has the same dimensions as r, r1, r2
rho_r = rho(ones(1,3),:,:);

r1_minus_robs = r1 - robs_r;
mag2_r1_minus_robs = sum(r1_minus_robs.^2,3);
mag_r1_minus_robs = sqrt(mag2_r1_minus_robs);

r2_minus_robs = r2 - robs_r;
mag2_r2_minus_robs = sum(r2_minus_robs.^2,3);
mag_r2_minus_robs = sqrt(mag2_r2_minus_robs);

r1_minus_rho = r1 - rho_r;
% mag_r1_minus_rho = sqrt(sum(r1_minus_rho.^2,3));

r2_minus_rho = r2 - rho_r;
% mag_r2_minus_rho = sqrt(sum(r2_minus_rho.^2,3));

%r1_minus_rho_dot_t = dot(r1_minus_rho, t, 3);
r1_minus_rho_dot_t = sum(r1_minus_rho.*t, 3);
%r2_minus_rho_dot_t = dot(r2_minus_rho, t, 3);
r2_minus_rho_dot_t = sum(r2_minus_rho.*t, 3);

%r1_minus_rho_dot_nu = dot(r1_minus_rho, nu, 3);
r1_minus_rho_dot_nu = sum(r1_minus_rho.*nu, 3);
%r2_minus_rho_dot_nu = dot(r2_minus_rho, nu, 3);
r2_minus_rho_dot_nu = sum(r2_minus_rho.*nu, 3);

abs_h = sqrt(h.*h); % abs() can not handle complex images
% abs_h matrix repeated - 3-by-N array, N is the triangle index
abs_h_r = abs_h(ones(1,3),:);

ln_arg_nom = mag_r2_minus_robs + r2_minus_rho_dot_t;
ln_arg_denom = mag_r1_minus_robs + r1_minus_rho_dot_t;

% If numerator and denominator are nearly zero (which is the case if r1-r2
% triangle edge and observation point lie on the same line) zero the
% logarithm value - it is multiplied by r1_minus_rho_dot_nu, which
% is zero in this case as well.
ln_arg_zeros = find(abs(ln_arg_nom) < 1e-15 | abs(ln_arg_denom) < 1e-15);
ln_arg_nom(ln_arg_zeros) = 1;
ln_arg_denom(ln_arg_zeros) = 1;

ln_arg = ln_arg_nom ./ ln_arg_denom;

atan1_arg = ( mag_r2_minus_robs + abs_h_r + r2_minus_rho_dot_t ) ...
               ./ r2_minus_rho_dot_nu;

atan2_arg = ( mag_r1_minus_robs + abs_h_r + r1_minus_rho_dot_t ) ...
               ./ r1_minus_rho_dot_nu;

% After summation over edges this gives one_over_r integral
one_over_r_edge = r1_minus_rho_dot_nu.*log(ln_arg) ...
                     - 2*abs_h_r.*(atan(atan1_arg) - atan(atan2_arg));

% If r1_minus_rho_dot_nu and r2_minus_rho_dot_nu are zero (if one is zero
% then another is zero also, this is the case when r1-r2 triangle edge
% and observation point projection to the triangle plane lie on the same
% line) then corresponding one_over_r_edge term is zero as well: logarithm
% is multiplied with r1_minus_rho_dot_nu, and both arctangents are pi because
% of infinite arguments. Enforce that to avoid problems caused by numerical
% jitter.
one_over_r_edge_zeros = find(abs(r1_minus_rho_dot_nu) < 1e-15);
one_over_r_edge(one_over_r_edge_zeros) = 0;

% Wrong sign in the paper?? It is corrected here
one_over_r = sum(one_over_r_edge, 1);

%r1_minus_robs_dot_t = dot(r1_minus_robs, t, 3);
r1_minus_robs_dot_t = sum(r1_minus_robs.*t, 3);
%r2_minus_robs_dot_t = dot(r2_minus_robs, t, 3);
r2_minus_robs_dot_t = sum(r2_minus_robs.*t, 3);

log1_arg = mag_r1_minus_robs + r1_minus_rho_dot_t;
log2_arg = mag_r2_minus_robs + r2_minus_rho_dot_t;

% Arguments can be zero if r1-r2 triangle edge and observation point
% lie on the same line, in this case the logarithm is multiplied by
% zero in the expression below, so it actual value does not matter.
log1_arg_zeros = find(abs(log1_arg) < 1e-15);
log1_arg(log1_arg_zeros) = 1;
log2_arg_zeros = find(abs(log2_arg) < 1e-15);
log2_arg(log2_arg_zeros) = 1;

v1 = mag_r2_minus_robs.*r2_minus_robs_dot_t ...
     - mag_r1_minus_robs.*r1_minus_robs_dot_t ...
      + ((mag2_r2_minus_robs - r2_minus_rho_dot_t.^2).*log(log2_arg)) ...
      - ((mag2_r1_minus_robs - r1_minus_rho_dot_t.^2).*log(log1_arg));
	   
v2 = 0.5*sum(nu.*v1(:,:,ones(1,3)), 1);

f_over_r = v2(ones(1,3),:,:) + (rho_r-r).*repmat(one_over_r, [3 1 3]);

% Cleanup intermediate variables to free up memory
clear robs_r r1_minus_robs mag2_r1_minus_robs mag2_r1_minus_robs;
clear r2_minus_robs mag2_r2_minus_robs mag_r2_minus_robs;
clear r1_minus_rho r2_minus_rho r1_minus_rho_dot_t r2_minus_rho_dot_t;
%clear r1_minus_rho_dot_nu ln_arg atan1_arg atan2_arg;
clear ln_arg_nom ln_arg_denom ln_arg_zeros;
clear r1_minus_robs_dot_t r2_minus_robs_dot_t v1 v2 one_over_r_edge;
clear one_over_r_edge_zeros;
clear log1_arg log1_arg_zeros log2_arg log2_arg_zeros;

% Number of quadrature points for the line integrals over edges.
qN=qn;

% Quadrature point positions and weights for the integrals over edges.
[qX,qW] = GLNodeWt(qN);
%[qX,qW] = GLTable(qN);

% Quadrature point positions. 3-by-N-by-3-by-qN array. First index is the
% local index of the egde in triangle, second one is the triangle index,
% third index is X-Y-Z, fourth is the quadrature point index.
qr = repmat(r1,[ 1 1 1 qN ]) + repmat(t,[ 1 1 1 qN ]) ...
        .* repmat(s,[ 1 1 3 qN ]) ...
          .* repmat(permute((qX*0.5 + 0.5),[ 2 3 4 1 ]), [ 3 N 3 1 ]);

% Distance from the observation point to quadrature point
qR = sqrt(sum((qr - repmat(robs,[ 3 1 1 qN ])).^2,3));

% Permute it so the third index is a quadrature point index.
qR = permute(qR,[1 2 4 3]);

% abs_h matrix repeated so it has the same dimensions as qR
abs_h_qr = abs_h_r(:,:,ones(1,qN));

% Temporary terms, exponent and sinc
ex = exp(-j*k*(qR+abs_h_qr)/2);
sinc_arg = k*(qR-abs_h_qr)/(2*pi);
ex_zeros = find(abs(ex) < 1e-15); % Is done to handle sinc overflow
sinc_arg(ex_zeros) = 1; % sinc of complex argument > +/-300i is Inf
si = sinc(sinc_arg);

% Calculate the quadrature over edges.
% Wrong sign in the paper?? It is corrected here
integrand1 = (ex.*si-1)./(qR+abs_h_qr);
qW_3 = repmat(permute(qW, [ 2 3 1 ]), [ 3 N ]);
integral1 = sum(integrand1.*qW_3,3).*(s/2);

e1_over_r = sum(r2_minus_rho_dot_nu.*integral1, 1);

e_over_r = e1_over_r + one_over_r;

% Cleanup intermediate variables to free up memory
%clear r2_minus_rho_dot_nu abs_h_qr ex si sinc_arg ex_zeros integral1;

% Vector from 'rho' point to the quadrature point.
q_rho = qr - repmat(rho,[ 3 1 1 qN ]);

% Squared distance from 'rho' point to the quadrature point
q_rho_2 = sum(q_rho.^2,3);

% Permute it so the third index is a quadrature point index.
q_rho_2 = permute(q_rho_2,[1 2 4 3]);

% Calculate the quadrature
integrand2 = integrand1.*q_rho_2;
integral2 = sum(integrand2.*qW_3,3).*(s/2);

v1 = nu.*repmat(integral2, [ 1 1 3 ]);
s1 = repmat(sum(v1,1), [ 3 1 1 ]);
f_e1_over_r = s1 + (rho_r - r).*repmat(e1_over_r, [ 3 1 3 ]);

f_e_over_r = f_e1_over_r + f_over_r;

% Cleanup intermediate variables to free up memory
%clear integrand1 q_rho q_rho_2 integrand2 integral2 v1 s1;

%===========================================================================
% f_cross_grad calculation begins here.
%
% f_cross_grad = f X grad(exp(-j*k*R)/R), where f is the RWG basis function
%  given by f = r - r(l), and r is the evaluation point coordinates and
%  r(l) is the free vertex.
%
% First, we compute the tangential (to the plane of triangle) part. Similar
% to the integrals above, we divide it into two parts:
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
integrand3 = exp(-j*k*qR)./qR;
integral3 = -h_r.*sum(integrand3.*qW_3,3).*(s/2);
fxg_obs = sum(nu.*repmat(integral3, [ 1 1 3 ]),1);
% Transform (x*dgdz,y*dgdz) into (y*dgdz,-x*dgdz)
fxg_obs = cross(fxg_obs, n, 3);

% Find indices of the elements for which h=0. Used later to zero the
% corresponding f_cross_grad elements.
%h_zeros = find(abs_h < 1e-15);

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
integrand4a = -exp(-j*k*qR)./(q_rho_2.*qR);
% The f_x and f_y given above is singular when rho=0 (and thus
% integrand4a is); to address that we add the following correction terms:
%   g_x = -x*z*exp(-j*k*abs(h))/(abs(h)*(x^2+y^2))
%   g_y = -y*z*exp(-j*k*abs(h))/(abs(h)*(x^2+y^2))
% so:
%   f_x_nonsing = f_x + g_x
%   f_y_nonsing = f_y + g_y
% Now it can be seen that both f_x_nonsing and f_y_nonsing are non-singular
% i.e. they converge to a finite limit when rho -> 0.
% This does not affect the integration results because divergence of vector
% function (g_x, g_y) is zero.
abs_h_qr_zeros = find(abs_h_qr < 1e-15);
abs_h_qr_fixed = abs_h_qr;
abs_h_qr_fixed(abs_h_qr_zeros) = 1e-15;
integrand4b = exp(-j*k*abs_h_qr_fixed)./(q_rho_2.*abs_h_qr_fixed);
integrand4 = integrand4a + integrand4b;
integral4 = h_r.*sum(integrand4.*qW_3,3).*(s/2);
fxg_r123 = r2_minus_rho_dot_nu.*integral4;
fxg_r123 = -(rho_r - r).*repmat(sum(fxg_r123,1), [ 3 1 3 ]);
% Transform (r(l)_x*dgdz,r(l)_y*dgdz) into (r(l)_y*dgdz,-r(l)_x*dgdz)
fxg_r123 = cross(fxg_r123, repmat(n, 3, 1), 3);

% The tangential part is done!
f_cross_grad_t = repmat(fxg_obs, 3, 1) - fxg_r123;

% Now compute the normal part. Similar to all the integrals above we divide
% it into two parts, the first part includes the observation coordinate only,
% and the second depends on the free vertex position.
% The first part is given by:
%  f1 = x*dgdy-y*dgdx
% f1 is always zero, so we only need to compute the integral of the second
% part which involves the free vertex position. The rest is fairly
% obvious I believe.
integrand5 = exp(-j*k*qR)./qR;
integral5 = sum(integrand5.*qW_3,3).*(s/2);
fxg_n = sum(nu.*repmat(integral5, [ 1 1 3 ]),1);
% Transform (dgdx,dgdy) into (dgdy,-dgdx)
fxg_n = cross(fxg_n, n, 3);
% This gives f(l)_x*dgdy-f(l)_y*dgdx
fxg_n = sum(repmat(fxg_n, 3, 1).*(rho_r - r), 3);

f_cross_grad_n = repmat(n,3,1).*repmat(fxg_n, [ 1 1 3 ]);

f_cross_grad = f_cross_grad_t + f_cross_grad_n;

% Diagnostics
badv = ~isfinite(f_e_over_r) | ~isfinite(f_cross_grad);
badidx = find(badv);
if ~isempty(badidx),
	error('i_calc : encountered inalid values');
end

