function test_tlfact
% Here we test if the factorized form of the transmission line Greens
% function is computed properly. We compute the factorized form coefficients
% and calculate the voltages/currents at a few points using these
% coefficients, after that we compute the came voltages/currents using
% "ordinary" Greens functions and compare the values obtained in both ways.

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


%
% Voltages and currests at all test points due to a current source
% at potint 32.
%
[ f0, f1, f2, f3, f4 ] = fact_vi(tl,3);

v31i_golden = calc_vi(tl,z31,3,z32,3);
v31i = calc_tlfact(tl, f0, f1, f2, f3, f4, 3, z31, z32);
assertEquals(v31i_golden, v31i, 1e-15);

v32i_golden = calc_vi(tl,z32,3,z32,3);
v32i = calc_tlfact(tl, f0, f1, f2, f3, f4, 3, z32, z32);
assertEquals(v32i_golden, v32i, 1e-15);

v33i_golden = calc_vi(tl,z33,3,z32,3);
v33i = calc_tlfact(tl, f0, f1, f2, f3, f4, 3, z33, z32);
assertEquals(v33i_golden, v33i, 1e-15);


[ f0, f1, f2, f3, f4 ] = fact_ii(tl, 3);

i31i_golden = calc_ii(tl,z31,3,z32,3);
% Notice the minus sign before f0, because in this case z<z'
i31i = calc_tlfact(tl, -f0, f1, f2, f3, f4, 3, z31, z32);
assertEquals(i31i_golden, i31i, 1e-7);

i33i_golden = calc_ii(tl,z33,3,z32,3);
i33i = calc_tlfact(tl, f0, f1, f2, f3, f4, 3, z33, z32);
assertEquals(i33i_golden, i33i, 1e-7);



%
% Voltages and currents at all test points due to a voltage source
% at potint 32.
%

[ f0, f1, f2, f3, f4 ] = fact_vv(tl,3);

v31v_golden = calc_vv(tl,z31,3,z32,3);
% Notice the minus sign before f0, because in this case z<z'
v31v = calc_tlfact(tl, -f0, f1, f2, f3, f4, 3, z31, z32);
assertEquals(v31v_golden, v31v, 1e-7);

v33v_golden = calc_vv(tl,z33,3,z32,3);
v33v = calc_tlfact(tl, f0, f1, f2, f3, f4, 3, z33, z32);
assertEquals(v33v_golden, v33v, 1e-7);


[ f0, f1, f2, f3, f4 ] = fact_iv(tl,3);

i31v_golden = calc_iv(tl,z31,3,z32,3);
i31v = calc_tlfact(tl, f0, f1, f2, f3, f4, 3, z31, z32);
assertEquals(i31v_golden, i31v, 1e-6);

i32v_golden = calc_iv(tl,z32,3,z32,3);
i32v = calc_tlfact(tl, f0, f1, f2, f3, f4, 3, z32, z32);
assertEquals(i32v_golden, i32v, 1e-6);

i33v_golden = calc_iv(tl,z33,3,z32,3);
i33v = calc_tlfact(tl, f0, f1, f2, f3, f4, 3, z33, z32);
assertEquals(i33v_golden, i33v, 1e-6);

