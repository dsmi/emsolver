function test_gspec
% Compare the spectral domain Greens functions calculated in two different
% ways - (1) using gspec_calc, which uses simplified formulas and does not
% use the transmission lines calculator (2) using gspec_validate, which
% uses the general case formulas and transmission lines calculator. Notice
% that gspec_validate does NOT work properly for kr (radial component of
% the wavevector) = 0
%


% Angular frequency
lay.freq = 1e10;

%
% First test - all three layers have different wavenumbers
%
lay.eps1 = eps0*5;
lay.eps2 = eps0*2;
lay.eps3 = eps0*7;
lay.z2 = 0.1;
lay.z3 = 0.2;

% Generate a number of kr/kz samples, to compare the G values at these
% points.
T02=5;
ns=100;
[ kr, kz ] = mkipath2(lay, T02, ns);

% gspec_validate doesn't work properly for kr=0 - remove the first sample
kr(1) = [];
kz(1,:) = [];

gs = gspec_calc(lay, kr, kz);
gspec_validate(gs, lay, kr, kz);

%
% Second test - layers 1 and 2 have the same wavenumbers, so
% the factors 2, 3 and 4 of the spectral domain Greens function
% should be zero.
%
lay.eps1 = eps0;
lay.eps2 = eps0;
lay.eps3 = eps0*7;
lay.z2 = 0.1;
lay.z3 = 0.2;

% Generate a number of kr/kz samples, to compare the G values at these
% points.
T02=5;
ns=100;
[ kr, kz ] = mkipath2(lay, T02,ns);

% gspec_validate doesn't work properly for kr=0 - remove the first sample
kr(1) = [];
kz(1,:) = [];

gs = gspec_calc(lay, kr, kz);
gspec_validate(gs, lay, kr, kz);

% Make sure factors 2, 3 and 4 are zero.
assertEquals(0, gs.Gvv2);
assertEquals(0, gs.Gvv3);
assertEquals(0, gs.Gvv4);

assertEquals(0, gs.Gzz2);
assertEquals(0, gs.Gzz3);
assertEquals(0, gs.Gzz4);

assertEquals(0, gs.Gzu2);
assertEquals(0, gs.Gzu3);
assertEquals(0, gs.Gzu4);

assertEquals(0, gs.Fvv2);
assertEquals(0, gs.Fvv3);
assertEquals(0, gs.Fvv4);

assertEquals(0, gs.Fzz2);
assertEquals(0, gs.Fzz3);
assertEquals(0, gs.Fzz4);

assertEquals(0, gs.Fzu2);
assertEquals(0, gs.Fzu3);
assertEquals(0, gs.Fzu4);

assertEquals(0, gs.Kf2);
assertEquals(0, gs.Kf3);
assertEquals(0, gs.Kf4);

assertEquals(0, gs.Cf2);
assertEquals(0, gs.Cf3);
assertEquals(0, gs.Cf4);

%
% Third test - layers 2 and 3 have the same wavenumbers, so
% the factors 1, 3 and 4 of the spectral domain Greens function
% should be zero.
%
lay.eps1 = eps0*7;
lay.eps2 = eps0;
lay.eps3 = eps0;
lay.z2 = 0.1;
lay.z3 = 0.2;

% Generate a number of kr/kz samples, to compare the G values at these
% points.
T02=5;
ns=100;
[ kr, kz ] = mkipath2(lay, T02, ns);

% gspec_validate doesn't work properly for kr=0 - remove the first sample
kr(1) = [];
kz(1,:) = [];

gs = gspec_calc(lay, kr, kz);
gspec_validate(gs, lay, kr, kz);

% Make sure factors 1, 3 and 4 are zero.
assertEquals(0, gs.Gvv1);
assertEquals(0, gs.Gvv3);
assertEquals(0, gs.Gvv4);

assertEquals(0, gs.Gzz1);
assertEquals(0, gs.Gzz3);
assertEquals(0, gs.Gzz4);

assertEquals(0, gs.Gzu1);
assertEquals(0, gs.Gzu3);
assertEquals(0, gs.Gzu4);

assertEquals(0, gs.Fvv1);
assertEquals(0, gs.Fvv3);
assertEquals(0, gs.Fvv4);

assertEquals(0, gs.Fzz1);
assertEquals(0, gs.Fzz3);
assertEquals(0, gs.Fzz4);

assertEquals(0, gs.Fzu1);
assertEquals(0, gs.Fzu3);
assertEquals(0, gs.Fzu4);

assertEquals(0, gs.Kf1);
assertEquals(0, gs.Kf3);
assertEquals(0, gs.Kf4);

assertEquals(0, gs.Cf1);
assertEquals(0, gs.Cf3);
assertEquals(0, gs.Cf4);

%
% Fourth test - layer 1 is a perfect conductor, layers 2 and 3 have the
% same wavenumbers.
%
lay.eps1 = Inf;
lay.eps2 = eps0;
lay.eps3 = eps0;
lay.z2 = 0.1;
lay.z3 = 0.2;

% Generate a number of kr/kz samples, to compare the G values at these
% points.
T02=5;
ns=100;
[ kr, kz ] = mkipath2(lay, T02, ns);

% gspec_validate doesn't work properly for kr=0 - remove the first sample
kr(1) = [];
kz(1,:) = [];

gs = gspec_calc(lay, kr, kz);
gspec_validate(gs, lay, kr,kz);

% Layers 2 and 3 have the same wavenumbers, so factors 1, 3 and 4
% should be zero.
assertEquals(0, gs.Gvv1);
assertEquals(0, gs.Gvv3);
assertEquals(0, gs.Gvv4);

assertEquals(0, gs.Gzz1);
assertEquals(0, gs.Gzz3);
assertEquals(0, gs.Gzz4);

assertEquals(0, gs.Gzu1);
assertEquals(0, gs.Gzu3);
assertEquals(0, gs.Gzu4);

assertEquals(0, gs.Fvv1);
assertEquals(0, gs.Fvv3);
assertEquals(0, gs.Fvv4);

assertEquals(0, gs.Fzz1);
assertEquals(0, gs.Fzz3);
assertEquals(0, gs.Fzz4);

assertEquals(0, gs.Fzu1);
assertEquals(0, gs.Fzu3);
assertEquals(0, gs.Fzu4);

assertEquals(0, gs.Kf1);
assertEquals(0, gs.Kf3);
assertEquals(0, gs.Kf4);

assertEquals(0, gs.Cf1);
assertEquals(0, gs.Cf3);
assertEquals(0, gs.Cf4);

%
% Fifth test - layer 3 is a perfect conductor, layers 1 and 2 have the
% same wavenumbers.
%
lay.eps1 = eps0;
lay.eps2 = eps0;
lay.eps3 = Inf;
lay.z2 = 0.1;
lay.z3 = 0.2;

% Generate a number of kr/kz samples, to compare the G values at these
% points.
T02=5;
ns=100;
[ kr, kz ] = mkipath2(lay, T02, ns);

% gspec_validate doesn't work properly for kr=0 - remove the first sample
kr(1) = [];
kz(1,:) = [];

gs = gspec_calc(lay, kr, kz);
gspec_validate(gs, lay, kr, kz);

% Layers 1 and 2 have the same wavenumbers, so factors 2, 3 and 4
% should be zero.
assertEquals(0, gs.Gvv2);
assertEquals(0, gs.Gvv3);
assertEquals(0, gs.Gvv4);

assertEquals(0, gs.Gzz2);
assertEquals(0, gs.Gzz3);
assertEquals(0, gs.Gzz4);

assertEquals(0, gs.Gzu2);
assertEquals(0, gs.Gzu3);
assertEquals(0, gs.Gzu4);

assertEquals(0, gs.Fvv2);
assertEquals(0, gs.Fvv3);
assertEquals(0, gs.Fvv4);

assertEquals(0, gs.Fzz2);
assertEquals(0, gs.Fzz3);
assertEquals(0, gs.Fzz4);

assertEquals(0, gs.Fzu2);
assertEquals(0, gs.Fzu3);
assertEquals(0, gs.Fzu4);

assertEquals(0, gs.Kf2);
assertEquals(0, gs.Kf3);
assertEquals(0, gs.Kf4);

assertEquals(0, gs.Cf2);
assertEquals(0, gs.Cf3);
assertEquals(0, gs.Cf4);
