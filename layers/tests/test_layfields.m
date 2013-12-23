function test_layfields
% 
% Compare fields of a dipole in layered media calculated in two different ways:
% (a) dyadic Greens function, Sommerfeld integrals evaluated numerically
% (b) mixed-potential representation, complex images
%
%

% Source dipole
l = [ 2e-7; 5e-7; 1e-7 ];

robs = [ -7.0e-3   2.0e-3   3.0e-3; ...  % r, observation position
          4.0e-3  -1.0e-3  -6.0e-3; ...
          1.5e-3   2.0e-3   2.5e-3 ];
rsrc = [  2.5e-3   1.0e-3  -5.0e-3; ...  % r', source position
         -3.2e-3   0.0e-3  -3.0e-3; ...
          2.5e-3   1.9e-3   1.5e-3 ];

% Several different stackups
lays(1).freq = 1e11;
lays(1).z3   = 3e-3;
lays(1).z2   = 1e-3;
lays(1).eps3 = eps0*2;
lays(1).eps2 = eps0*10;
lays(1).eps1 = eps0*7;
lays(2).freq = 1e11;
lays(2).z3   = 3e-3;
lays(2).z2   = 1e-3;
lays(2).eps3 = eps0*3;
lays(2).eps2 = eps0*10;
lays(2).eps1 = eps0*5;
lays(3).freq = 1e11;
lays(3).z3   = 3e-3;
lays(3).z2   = 1e-3;
lays(3).eps3 = eps0;
lays(3).eps2 = eps0*10;
lays(3).eps1 = eps0-j*Inf;
lays(4).freq = 1e11;
lays(4).z3   = 3e-3;
lays(4).z2   = 1e-3;
lays(4).eps3 = eps0-j*Inf;
lays(4).eps2 = eps0*10;
lays(4).eps1 = eps0*4;

for il=1:length(lays),
	lay = lays(il);
	gi = mkimages(lay);
	for ir=1:size(robs,2),
		ro = robs(:,ir);
		rs = rsrc(:,ir);
		[ efj1 efm1 ] = dfield_dyg(lay,ro,2,rs,2,l);
		[ efj2 efm2 ] = dfield_dcim(lay,gi,ro,rs,l);
		jreltol(il,ir) = norm(efj1-efj2)/norm(efj1);
		mreltol(il,ir) = norm(efm1-efm2)/norm(efm1);
	end
end

assertTrue(~nnz(jreltol>0.02))
assertTrue(~nnz(mreltol>0.02))

