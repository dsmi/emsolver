function test_tlines
% Test if the tlines calculator computes the voltages and currents
% properly. Golden values have been obtained using the eldo simulator.
% The test fixture consists of five various tlines, open at ends.
% The voltage/current source is in the 3rd (middle) tline, observation
% points are at the left and right in the same tline and one observation
% point in each of the other tlines.

% Tline parameters
R   = [ 10    5      4     8     3     ]; 
L   = [ 1e-6  1.5e-6 8e-7  5e-7  1e-7  ];
G   = [ 5e2   4e2    5e1   1e2   2e2   ];
C   = [ 1e-11 3e-11  5e-11 5e-12 4e-11 ];
len = [ 1e-3  4e-3   4e-4  2e-3  6e-3  ];

% Frequency and angular frequency
freq=1e6;
afreq=2*pi*freq;

% Characteristic impedances
Z0=sqrt( (R+j*afreq*L)./(G+j*afreq*C) );

% Propagation constants
tl_k=sqrt( (R+j*afreq*L).*(G+j*afreq*C) );

% Output the tline parameters
if 0,
	for l=1:length(R),
		printf("Tline %i\n",l);
		printf("    R=%e, L=%e, G=%e, c=%e\n", R(l), L(l), G(l), C(l));
		printf("    len=%e\n", len(l));
		printf("    k=%e + %ei\n", real(tl_k(l)), imag(tl_k(l)));
		printf("    Z0=%e + %ei\n", real(Z0(l)), imag(Z0(l)));
	end
	printf("\n");
end

% Output the matched terminations
if 0,
	% Matched terminations
	Rterm1=real(Z0(1));
	Lterm1=imag(Z0(1))/afreq;
	Rterm5=real(Z0(5));
	Lterm5=imag(Z0(5))/afreq;

	printf("Matched termination for line 1: R=%e, L=%e\n", Rterm1, Lterm1);
	printf("Matched termination for line 5: R=%e, L=%e\n", Rterm5, Lterm5);
	printf("\n");
end

% Setup the endpoint coordinates to be used by calculator.
tl_z = zeros(1,length(len)+1);
for n=2:length(len)+1,
	tl_z(n)=sum(len(1:n-1));
end

% Matched termination at both ends
Gls1=0;
GgrN=0;

% Run the coefficients precomputation
tl = calc_tlines(tl_z, Z0, tl_k, Gls1, GgrN);

% Coordinate of the observation point in the tline1
z11 = tl.z(1)+len(1)*0.2;

% Coordinate of the observation point in the tline2
z21 = tl.z(2)+len(2)*0.7;

% Coordinate of the observation point in the source tline at the left
% of the source.
z31 = tl.z(3)+len(3)*0.1;

% Coordinate of the source.
z32 = tl.z(3)+len(3)*0.4;

% Coordinate of the observation point in the source tline at the right
% of the source.
z33 = tl.z(3)+len(3)*0.8;

% Coordinate of the observation point in the tline4
z41 = tl.z(4)+len(4)*0.8;

% Coordinate of the observation point in the tline5
z51 = tl.z(5)+len(2)*0.3;

if 0,
	% To insert the source and the probe points in the netlist, we split
	% tlines into parts. There are the part lenghts.
	len1a=z11-tl.z(1)
	len1b=len(1)-len1a

	len2a=z21-tl.z(2)
	len2b=len(2)-len2a

	len3a=z31-tl.z(3)
	len3b=z32-z31
	len3c=z33-z32
	len3d=tl.z(4)-z33

	len4a=z41-tl.z(4)
	len4b=len(4)-len4a

	len5a=z51-tl.z(5)
	len5b=len(5)-len5a
end

%
% Voltages and currests at all test points due to a current source
% at potint 32.
%
v11i_golden = 5.530704670E-002 + j*2.792725160E-003;
v11i = calc_vi(tl,z11,1,z32,3);
assertEquals(v11i_golden, v11i, 1e-8);

v21i_golden = 6.635229415E-002 + j*1.364712873E-002;
v21i = calc_vi(tl,z21,2,z32,3);
assertEquals(v21i_golden, v21i, 1e-8);

v31i_golden = 6.996286751E-002 + j*1.852721015E-002;
v31i = calc_vi(tl,z31,3,z32,3);
assertEquals(v31i_golden, v31i, 1e-7);

v32i_golden = 7.022983512E-002 + j*1.878425335E-002;
v32i = calc_vi(tl,z32,3,z32,3);
assertEquals(v32i_golden, v32i, 1e-7);

v33i_golden = 6.994600108E-002 + j*1.832320836E-002;
v33i = calc_vi(tl,z33,3,z32,3);
assertEquals(v33i_golden, v33i, 1e-7);

v41i_golden = 6.349720628E-002 + j*1.470421942E-002;
v41i = calc_vi(tl,z41,4,z32,3);
assertEquals(v41i_golden, v41i, 1e-7);

v51i_golden = 6.017657551E-002 + j*1.328212361E-002;
v51i = calc_vi(tl,z51,5,z32,3);
assertEquals(v51i_golden, v51i, 1e-7);


i11i_golden = -3.508310054E-001 + j*8.215917150E-002;
i11i = calc_ii(tl,z11,1,z32,3);
assertEquals(i11i_golden, i11i, 1e-7);

i21i_golden = -4.435226740E-001 + j*7.113700823E-002;
i21i = calc_ii(tl,z21,2,z32,3);
assertEquals(i11i_golden, i11i, 1e-7);

i31i_golden = -4.763524250E-001 + j*6.341536060E-002;
i31i = calc_ii(tl,z31,3,z32,3);
assertEquals(i31i_golden, i31i, 1e-7);

i33i_golden = 5.226662954E-001 + j*6.315499071E-002;
i33i = calc_ii(tl,z33,3,z32,3);
assertEquals(i33i_golden, i33i, 1e-7);

i41i_golden = 5.117243298E-001 + j*6.045954683E-002;
i41i = calc_ii(tl,z41,4,z32,3);
assertEquals(i41i_golden, i41i, 1e-7);

i51i_golden = 4.945627246E-001 + j*5.663023209E-002;
i51i = calc_ii(tl,z51,5,z32,3);
assertEquals(i51i_golden, i51i, 1e-7);


%
% Voltages and currents at all test points due to a voltage source
% at potint 32.
%
v11v_golden = -3.9982741858E-001 + j*3.6282370763E-002;
v11v = calc_vv(tl,z11,1,z32,3);
assertEquals(v11v_golden, v11v, 1e-7);

v21v_golden = -4.9016281360E-001 - j*3.0379269743E-002;
v21v = calc_vv(tl,z21,2,z32,3);
assertEquals(v21v_golden, v21v, 1e-7);

v31v_golden = -5.210490E-001 - j*6.173031E-002;
v31v = calc_vv(tl,z31,3,z32,3);
assertEquals(v31v_golden, v31v, 1e-7);

v33v_golden = 4.738695E-001 - j*6.540092E-002;
v33v = calc_vv(tl,z33,3,z32,3);
assertEquals(v33v_golden, v33v, 1e-7);

v41v_golden = 4.2528720098E-001 - j*7.1162323481E-002;
v41v = calc_vv(tl,z41,4,z32,3);
assertEquals(v41v_golden, v41v, 1e-7);

v51v_golden = 4.0139034210E-001 - j*7.1431795034E-002;
v51v = calc_vv(tl,z51,5,z32,3);
assertEquals(v51v_golden, v51v, 1e-7);


i11v_golden = 2.4345203745E+000 - j*9.4702970954E-001;
i11v = calc_iv(tl,z11,1,z32,3);
assertEquals(i11v_golden, i11v, 1e-6);

i21v_golden = 3.1110686329E+000 - j*9.6231752171E-001;
i21v = calc_iv(tl,z21,2,z32,3);
assertEquals(i21v_golden, i21v, 1e-6);

i31v_golden = 3.354578E+000 - j*9.403288E-001;
i31v = calc_iv(tl,z31,3,z32,3);
assertEquals(i31v_golden, i31v, 1e-6);

i32v_golden = 3.3577111637E+000 - j*9.3995371091E-001;
i32v = calc_iv(tl,z32,3,z32,3);
assertEquals(i32v_golden, i32v, 1e-6);

i33v_golden = 3.3539085922E+000 - j*9.3943891562E-001;
i33v = calc_iv(tl,z33,3,z32,3);
assertEquals(i33v_golden, i33v, 1e-6);

i41v_golden = 3.2802125169E+000 - j*9.2816344049E-001;
i41v = calc_iv(tl,z41,4,z32,3);
assertEquals(i41v_golden, i41v, 1e-6);

i51v_golden = 3.1656359923E+000 - j*9.0804508281E-001;
i51v = calc_iv(tl,z51,5,z32,3);
assertEquals(i51v_golden, i51v, 1e-6);
