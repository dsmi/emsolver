function test_roundz
%
% "Overclocked capacitor" - impedance of the waveguide formed by a pair of
% round plates. Compare the results against the analytical solution.
%

% The geometry - round with a hole in center which is the port.
r1 = 0.001;
r2 = 1;
d = 0.01;
n1 = 10;
n2 = 40;

[ e1, v1 ] = mkcir2d(r1, n1);
e1(:, [ 1 2 ]) = e1(:, [ 2 1 ]);

[ e2, v2 ] = mkcir2d(r2, n2);
e = [ e1; e2 + size(v1, 1) ];
v = [ v1; v2 ];

% The port edges
c1 = find_edges2d(e, v, 0, 0, r1*1.1);
ports = { c1' };

%plotmesh2d(e,v,ports,0);

freqs = logspace(6,9,20);
Za = []; % Analytical solution
Z1 = []; % Simulated
Z2 = []; % Simulated via admittance

for freq = freqs,

	% Parameters of the plane.
	Yplane = j*freq*eps0/d;
	Zplane = j*freq*mu0*d;
	k = sqrt(-Yplane*Zplane);

	fintgsl = @(rsrc,robs) intg_helmsl2d(k,rsrc,robs);
	fintgdl = @(rsrc,robs) intg_helmdl2d(k,rsrc,robs);

	Z = extractz2(e, v, ports, fintgsl, fintgdl)*Zplane;
	Z1 = [ Z1 Z ];

	Y=extracty2(e, v, ports, fintgsl, fintgdl)/Zplane;
	Z=1/Y;
	Z2 = [ Z2 Z ];

	%C = eps0*pi*r2*r2/d;
	%Zc = 1/(j*freq*C);

	% Analytical solution
	M = [ besselj(0, k*r1) bessely(0, k*r1) ; ...
		  besselj(0, k*r2) bessely(0, k*r2) ];
	%M = [ besselh(0,2,k*r1) besselh(0,1,k*r1) ; ...
	%      besselh(0,2,k*r2) besselh(0,1,k*r2) ]

	dM = [ -besselj(1, k*r1)*k -bessely(1, k*r1)*k ; ...
		   -besselj(1, k*r2)*k -bessely(1, k*r2)*k ];
	%dM = [ -besselh(1,2,k*r1)*k -besselh(1,1,k*r1)*k ; ...
	%       -besselh(1,2,k*r2)*k -besselh(1,1,k*r2)*k ]

	b = [ 1 ; 0 ];
	x = dM\b;

	V = M*x;
	I = -2*pi*r1/Zplane;
	Z = V(1)/I;
	Za = [ Za Z ];
end

assertEquals(Z1, Z2, Z1*1e-10);
assertEquals(Za, Z1, Za*3e-2);


