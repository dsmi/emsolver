function test_manytl
% Test if the tlines calculator works properly if given a number of
% chains of tlines at once. 

% Tline parameters
R1   = [ 10    5      4     8     3     ]; 
L1   = [ 1e-6  1.5e-6 8e-7  5e-7  1e-7  ];
G1   = [ 5e2   4e2    5e1   1e2   2e2   ];
C1   = [ 1e-11 3e-11  5e-11 5e-12 4e-11 ];
len1 = [ 1e-3  4e-3   4e-4  2e-3  6e-3  ];
R_all   = [ R1   ; R1*0.2 ; R1*1.1 ]; 
L_all   = [ L1   ; L1*2.0 ; L1*0.9 ];
G_all   = [ G1   ; G1*0.8 ; G1*0.6 ];
C_all   = [ C1   ; G1*0.5 ; G1*1.3 ];
len_all = [ len1 ; len1   ; len1   ];

% Frequency and angular frequency
freq=1e6;
afreq=2*pi*freq;

% Characteristic impedances
Z0_all=sqrt( (R_all+j*afreq*L_all)./(G_all+j*afreq*C_all) );

% Propagation constants
k_all=sqrt( (R_all+j*afreq*L_all).*(G_all+j*afreq*C_all) );

% Tline endpoint coordinates
z_all = zeros(3,length(len1)+1);
for n=2:length(len1)+1,
	z_all(:,n)=sum(len_all(:,1:n-1),2);
end

% Termination paramaters
Gls1_all = [ 1 ; -1 ;  0 ];
GgrN_all = [ 0 ;  1 ; -1 ];

% Coordinates and layers of the observation points
obsz = [ z_all(1,1)+len1(1)*0.2; z_all(1,2)+len1(2)*0.7; z_all(1,3)+len1(3)*0.1; ...
         z_all(1,3)+len1(3)*0.8; z_all(1,4)+len1(4)*0.8; z_all(1,5)+len1(2)*0.3 ];
obsj = [ 1 2 3 3 4 5 ];

% Coordinate and layer of the source.
zsrc = z_all(1,3)+len1(3)*0.4;
jsrc = 3;

% Compute the test values by processing each chain individually
test_vi = zeros(3,length(obsz));
test_ii = zeros(3,length(obsz));
test_vv = zeros(3,length(obsz));
test_iv = zeros(3,length(obsz));

for tl_idx=1:3,
	% Setup parameters for the current chain
	tl_z = z_all(tl_idx,:);
	Z0 = Z0_all(tl_idx,:);
	tl_k = k_all(tl_idx,:);
	Gls1 = Gls1_all(tl_idx);
	GgrN = GgrN_all(tl_idx);
	% Coefficients precomputation
	tl=calc_tlines(tl_z, Z0, tl_k, Gls1, GgrN);
	for obs_idx=1:length(obsz),
		zobs = obsz(obs_idx);
		jobs = obsj(obs_idx);
		test_vi(tl_idx,obs_idx) = calc_vi(tl,zobs,jobs,zsrc,jsrc);
		test_ii(tl_idx,obs_idx) = calc_ii(tl,zobs,jobs,zsrc,jsrc);
		test_vv(tl_idx,obs_idx) = calc_vv(tl,zobs,jobs,zsrc,jsrc);
		test_iv(tl_idx,obs_idx) = calc_iv(tl,zobs,jobs,zsrc,jsrc);
	end
end

% Setup parameters for all chains at once
tl_z = z_all;
Z0 = Z0_all;
tl_k = k_all;
Gls1 = Gls1_all;
GgrN = GgrN_all;

% Coefficients precomputation
tl = calc_tlines(tl_z, Z0, tl_k, Gls1, GgrN);

vi = [];
ii = [];
vv = [];
iv = [];
for obs_idx=1:length(obsz),
%for obs_idx=3:3,
	zobs = obsz(obs_idx);
	jobs = obsj(obs_idx);
	vi = [ vi calc_vi(tl,zobs,jobs,zsrc,jsrc) ];
	ii = [ ii calc_ii(tl,zobs,jobs,zsrc,jsrc) ];
	vv = [ vv calc_vv(tl,zobs,jobs,zsrc,jsrc) ];
	iv = [ iv calc_iv(tl,zobs,jobs,zsrc,jsrc) ];
end

assertEquals(test_vi,vi);
assertEquals(test_ii,ii);
assertEquals(test_vv,vv);
assertEquals(test_iv,iv);
