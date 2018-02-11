function test_ilay
% 
% Test the ilay_fp and ilay_p functions by comparison with quadrature.
% 

% The stackup
lay.freq = 1e11;
lay.z3   = 5e-4;
lay.z2   = -5e-4;
lay.eps3 = eps0*3;
lay.eps2 = eps0*10;
lay.eps1 = eps0*7;

gi = mkimages(lay);

% Testing primitive - a box.
[ tri, x, y, z ] = mkbox(1e-3, 1.2e-3, 6e-4, 1, 1, 1);

% Rotate it a little
[ x, y, z ] = rotmesh(x, y, z, 0.1, 0.3, -0.2);

% Number of observation points/source triangles
N = size(tri, 1);

% The source triangles
sx = x(tri);
sy = y(tri);
sz = z(tri);

% The observation points
ox = ((1:N)*3e-3/N)';
oy = -((1:N)*4e-4/N)';
oz = ((1:N)*4e-4/N)';

rsrc = cat(3, sx, sy, sz);
robs = [ ox oy oz ];

[ fp_test p_test ] = ilay_q(lay, gi, rsrc, robs, 64);

m = 1/(-4*pi*j*lay.freq*eps0);
p = ilay_p(lay, gi, rsrc, robs, 64)*m;

m = j*lay.freq*mu0/(4*pi); % This multiplier is applied to the integration
                           %  result of ilay_fp in mklhsmat
fp = ilay_fp(lay, gi, rsrc, robs, 64)*m;

assertTrue(~nnz(p_test-p > 1e-15))
assertTrue(~nnz(fp_test-fp > 2e-15))

m = j*lay.eps2*lay.freq/(4*pi);
fxg = ilay_fxg(lay, gi, rsrc, robs, 64)*m;
fxg_test = ilay_q2(lay, gi, rsrc, robs, 64);

assertTrue(~nnz(fxg-fxg_test > 1e-15))
