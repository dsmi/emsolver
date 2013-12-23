function test_reflect
% Test if the tlines calculator computes the reflection and transfer
% coefficients properly. Golden values have been obtained using the 
% eldo simulator. The test fixture consists of five various tlines,
% matched, shorted or open at ends.

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

% The reference values of the admittances (please mention that these are full
% admittances which include both forward and backward traveling waves and
% are different from characteristic admittances) have been obtained using
% eldo simulator, see reflect1.cir in netlists directory.
Y2=6.887481E+000-i*2.317741E+000;
Y3=7.141773E+000-i*1.049352E+000;
Y4=7.231980E+000-i*9.707909E-001;
Y5=8.034806E+000-i*8.323680E-001;
Y6=8.034807E+000-i*8.323680E-001;

% Calculate golden reflection coefficients based on the reference admittances
% calculated by eldo.
tl.Ggr_1_golden=(tl.Y0(1)-Y2)/(tl.Y0(1)+Y2);
tl.Ggr_2_golden=(tl.Y0(2)-Y3)/(tl.Y0(2)+Y3);
tl.Ggr_3_golden=(tl.Y0(3)-Y4)/(tl.Y0(3)+Y4);
tl.Ggr_4_golden=(tl.Y0(4)-Y5)/(tl.Y0(4)+Y5);
tl.Ggr_5_golden=(tl.Y0(5)-Y6)/(tl.Y0(5)+Y6);

% Compare the calculated reflection coefficients with the golden.
assertEquals(tl.Ggr_1_golden, tl.Ggr(1), 1e-6);
assertEquals(tl.Ggr_2_golden, tl.Ggr(2), 1e-6);
assertEquals(tl.Ggr_3_golden, tl.Ggr(3), 1e-6);
assertEquals(tl.Ggr_4_golden, tl.Ggr(4), 1e-6);
assertEquals(tl.Ggr_5_golden, tl.Ggr(5), 1e-6);

% The reference values of the admittances (please mention that these are full
% admittances which include both forward and backward traveling waves and
% are different from characteristic admittances) have been obtained using
% eldo simulator, see reflect2.cir in netlists directory.
Y1=-6.252380E+000+i*1.801223E+000;
Y2=-6.252380E+000+i*1.801223E+000;
Y3=-6.147429E+000+i*2.530569E+000;
Y4=-6.055916E+000+i*2.542949E+000;
Y5=-5.627309E+000+i*2.271797E+000;

% Calculate golden reflection coefficients based on the reference admittances
% calculated by eldo.
tl.Gls_1_golden=(tl.Y0(1)+Y1)/(tl.Y0(1)-Y1);
tl.Gls_2_golden=(tl.Y0(2)+Y2)/(tl.Y0(2)-Y2);
tl.Gls_3_golden=(tl.Y0(3)+Y3)/(tl.Y0(3)-Y3);
tl.Gls_4_golden=(tl.Y0(4)+Y4)/(tl.Y0(4)-Y4);
tl.Gls_5_golden=(tl.Y0(5)+Y5)/(tl.Y0(5)-Y5);

% Compare the calculated reflection coefficients with the golden.
assertEquals(tl.Gls_1_golden, tl.Gls(1), 1e-6);
assertEquals(tl.Gls_2_golden, tl.Gls(2), 1e-6);
assertEquals(tl.Gls_3_golden, tl.Gls(3), 1e-6);
assertEquals(tl.Gls_4_golden, tl.Gls(4), 1e-6);
assertEquals(tl.Gls_5_golden, tl.Gls(5), 1e-6);

% The reference transfer coefficients, also from eldo simulation.
tl.Tls1_golden=9.286090E-001-i*1.975693E-002;
tl.Tls2_golden=7.982089E-001-i*1.539854E-001;
tl.Tls3_golden=9.852136E-001-i*8.087461E-003;
tl.Tls4_golden=8.972334E-001+i*1.598661E-003;
tl.Tls5_golden=8.917780E-001+i*1.383757E-002;

% Validate the calculated transfer coefficients.
assertEquals(tl.Tls1_golden, tl.Tls(1), 1e-6);
assertEquals(tl.Tls2_golden, tl.Tls(2), 1e-6);
assertEquals(tl.Tls3_golden, tl.Tls(3), 1e-6);
assertEquals(tl.Tls4_golden, tl.Tls(4), 1e-6);
assertEquals(tl.Tls5_golden, tl.Tls(5), 1e-6);


% Next, short termination at both ends
Gls1=-1;
GgrN=-1;

% Run the coefficients precomputation
tl = calc_tlines(tl_z, Z0, tl_k, Gls1, GgrN);

% The reference values of the admittances (please mention that these are full
% admittances which include both forward and backward traveling waves and
% are different from characteristic admittances) have been obtained using
% eldo simulator, see reflect3.cir in netlists directory.
Y2=1.087972E+001-i*9.191697E+000;
Y3=2.541252E+001-i*8.575276E+000;
Y4=2.727133E+001-i*8.020454E+000;
Y5=5.362045E+001-i*1.114670E+001;

% Calculate golden reflection coefficients based on the reference admittances
% calculated by eldo.
tl.Ggr_1_golden=(tl.Y0(1)-Y2)/(tl.Y0(1)+Y2);
tl.Ggr_2_golden=(tl.Y0(2)-Y3)/(tl.Y0(2)+Y3);
tl.Ggr_3_golden=(tl.Y0(3)-Y4)/(tl.Y0(3)+Y4);
tl.Ggr_4_golden=(tl.Y0(4)-Y5)/(tl.Y0(4)+Y5);
tl.Ggr_5_golden=-1;

% Compare the calculated reflection coefficients with the golden.
assertEquals(tl.Ggr_1_golden, tl.Ggr(1), 1e-6);
assertEquals(tl.Ggr_2_golden, tl.Ggr(2), 1e-6);
assertEquals(tl.Ggr_3_golden, tl.Ggr(3), 1e-6);
assertEquals(tl.Ggr_4_golden, tl.Ggr(4), 1e-6);
assertEquals(tl.Ggr_5_golden, tl.Ggr(5), 1e-6);

% The reference values of the admittances (please mention that these are full
% admittances which include both forward and backward traveling waves and
% are different from characteristic admittances) have been obtained using
% eldo simulator, see reflect4.cir in netlists directory.
Y2=-7.186229E+001+i*4.504776E+001;
Y3=-1.124580E+001+i*1.559292E+001;
Y4=-1.076910E+001+i*1.483568E+001;
Y5=-1.006132E+001+i*1.037728E+001;

% Calculate golden reflection coefficients based on the reference admittances
% calculated by eldo.
tl.Gls_1_golden=-1;
tl.Gls_2_golden=(tl.Y0(2)+Y2)/(tl.Y0(2)-Y2);
tl.Gls_3_golden=(tl.Y0(3)+Y3)/(tl.Y0(3)-Y3);
tl.Gls_4_golden=(tl.Y0(4)+Y4)/(tl.Y0(4)-Y4);
tl.Gls_5_golden=(tl.Y0(5)+Y5)/(tl.Y0(5)-Y5);

% Compare the calculated reflection coefficients with the golden.
assertEquals(tl.Gls_1_golden, tl.Gls(1), 1e-6);
assertEquals(tl.Gls_2_golden, tl.Gls(2), 1e-6);
assertEquals(tl.Gls_3_golden, tl.Gls(3), 1e-6);
assertEquals(tl.Gls_4_golden, tl.Gls(4), 1e-6);
assertEquals(tl.Gls_5_golden, tl.Gls(5), 1e-6);

% The reference transfer coefficients, also from eldo simulation.
tl.Tls1_golden=0.0;
tl.Tls2_golden=1.999490E-001-i*9.057060E-002;
tl.Tls3_golden=9.529563E-001+i*2.104333E-003;
tl.Tls4_golden=7.752746E-001+i*1.034555E-001;
tl.Tls5_golden=8.003322E-001+i*9.547254E-002;

% Validate the calculated transfer coefficients.
assertEquals(tl.Tls1_golden, tl.Tls(1), 1e-7);
assertEquals(tl.Tls2_golden, tl.Tls(2), 1e-7);
assertEquals(tl.Tls3_golden, tl.Tls(3), 1e-7);
assertEquals(tl.Tls4_golden, tl.Tls(4), 1e-7);
assertEquals(tl.Tls5_golden, tl.Tls(5), 1e-7);


% Finally, open termination at both ends
Gls1=1;
GgrN=1;

% Run the coefficients precomputation
tl = calc_tlines(tl_z, Z0, tl_k, Gls1, GgrN);

% The reference values of the admittances (please mention that these are full
% admittances which include both forward and backward traveling waves and
% are different from characteristic admittances) have been obtained using
% eldo simulator, see reflect3.cir in netlists directory.
Y2=2.871867E+000-i*1.869732E-001;
Y3=1.382021E+000-i*1.552758E-002;
Y4=1.365115E+000-i*1.179461E-002;
Y5=1.191431E+000-i*1.777201E-003;

% Calculate golden reflection coefficients based on the reference admittances
% calculated by eldo.
tl.Ggr_1_golden=(tl.Y0(1)-Y2)/(tl.Y0(1)+Y2);
tl.Ggr_2_golden=(tl.Y0(2)-Y3)/(tl.Y0(2)+Y3);
tl.Ggr_3_golden=(tl.Y0(3)-Y4)/(tl.Y0(3)+Y4);
tl.Ggr_4_golden=(tl.Y0(4)-Y5)/(tl.Y0(4)+Y5);
tl.Ggr_5_golden=1;

% Compare the calculated reflection coefficients with the golden.
assertEquals(tl.Ggr_1_golden, tl.Ggr(1), 1e-7);
assertEquals(tl.Ggr_2_golden, tl.Ggr(2), 1e-7);
assertEquals(tl.Ggr_3_golden, tl.Ggr(3), 1e-7);
assertEquals(tl.Ggr_4_golden, tl.Ggr(4), 1e-7);
assertEquals(tl.Ggr_5_golden, tl.Ggr(5), 1e-7);

% The reference values of the admittances (please mention that these are full
% admittances which include both forward and backward traveling waves and
% are different from characteristic admittances) have been obtained using
% eldo simulator, see reflect4.cir in netlists directory.
Y2=-4.991677E-001+i*5.214473E-004;
Y3=-2.059169E+000+i*6.905508E-002;
Y4=-2.071745E+000+i*7.713873E-002;
Y5=-2.196552E+000+i*9.969568E-002;

% Calculate golden reflection coefficients based on the reference admittances
% calculated by eldo.
tl.Gls_1_golden=1;
tl.Gls_2_golden=(tl.Y0(2)+Y2)/(tl.Y0(2)-Y2);
tl.Gls_3_golden=(tl.Y0(3)+Y3)/(tl.Y0(3)-Y3);
tl.Gls_4_golden=(tl.Y0(4)+Y4)/(tl.Y0(4)-Y4);
tl.Gls_5_golden=(tl.Y0(5)+Y5)/(tl.Y0(5)-Y5);

% Compare the calculated reflection coefficients with the golden.
assertEquals(tl.Gls_1_golden, tl.Gls(1), 1e-7);
assertEquals(tl.Gls_2_golden, tl.Gls(2), 1e-7);
assertEquals(tl.Gls_3_golden, tl.Gls(3), 1e-7);
assertEquals(tl.Gls_4_golden, tl.Gls(4), 1e-7);
assertEquals(tl.Gls_5_golden, tl.Gls(5), 1e-7);

% The reference transfer coefficients, also from eldo simulation.
tl.Tls1_golden=9.975032E-001-i*1.564269E-003;
tl.Tls2_golden=9.726403E-001-i*4.677424E-002;
tl.Tls3_golden=9.965461E-001-i*4.022007E-003;
tl.Tls4_golden=9.658126E-001-i*1.159156E-002;
tl.Tls5_golden=9.515244E-001-i*7.977060E-003;

% Validate the calculated transfer coefficients.
assertEquals(tl.Tls1_golden, tl.Tls(1), 1e-7);
assertEquals(tl.Tls2_golden, tl.Tls(2), 1e-7);
assertEquals(tl.Tls3_golden, tl.Tls(3), 1e-7);
assertEquals(tl.Tls4_golden, tl.Tls(4), 1e-7);
assertEquals(tl.Tls5_golden, tl.Tls(5), 1e-7);
