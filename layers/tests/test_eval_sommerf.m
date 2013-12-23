function test_eval_sommerf
%
% Test the eval_sommerf function which evaluates the Sommerfeld integral
% numerically. 
% In this test the spectral domain representation of the free-space Green's
% function is transformed to the spatial domain in two ways: using
% eval_sommerf function and using the Sommerfeld identity. In other words,
% the result obtained via the Sommerfeld identity is compared with
% the result obtained using the numerical integration.
%  The Sommerfeld identity is:
%     exp(-i*k*R)/R = \int_0^inf exp(-i*kz*|z-z'|)/i*kz*J0(kr*rho)*kr*dkr
%  Where:
%        R = |r-r'|
%        J0 is the Bessel function
%

% Angular frequency
lay.freq = 1e10;

% Setup the layers stackup
lay.z2   = 1e-3;
lay.z3   = 3e-3;
lay.eps1 = eps0;
lay.eps2 = eps0;
lay.eps3 = eps0;

dx = 5e-2;
dy = 7e-2;
dz = 9e-2;

r = sqrt(dx*dx+dy*dy+dz*dz);
rho = sqrt(dx*dx+dy*dy);

lay_k = calc_lay_k(lay);
ki = lay_k(2);

fg = @(kr, kz) exp(-j.*kz(:,2)*abs(dz))./(j*kz(:,2)).*kr;
somm = eval_sommerf(lay, 0, fg, rho);

% Test value. eval_sommerf applies 1/(2*pi) factor
testv = 1/(2*pi)*exp(-j*ki*r)/r;

assertEquals(testv, somm, 1e-5);
