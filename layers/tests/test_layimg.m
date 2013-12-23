function test_layimg
% 
% Test if the Greens function for the layered media give correct fields in
% the case of perfectly conducting plane on top or on bottom. In this case
% we can calculate the fields using the image theory and compare the fields
% calculated using the layered Greens function+dcim against it.
%

% Angular frequency
freq = lay.freq = 1e11;

% Source dipole
l = [ 2e-7; 5e-7; 1e-7 ];

robs = [ -7.0e-3; 4.0e-3; 1.5e-3 ]; % r, observation position
rsrc = [ 2.5e-3; -3.2e-3; 2.5e-3 ]; % r', source position

% First case - perfectly conducting plane on bottom
lay.z3   = 3e-3;
lay.z2   = 1e-3;
lay.eps3 = eps0;
lay.eps2 = eps0;
lay.eps1 = eps0-j*Inf;

gi = mkimages(lay);

ef1 = dfield_dcim(lay,gi,robs,rsrc,l);

[ dx, dy, dz ] = uncat(1, rsrc-robs);
ll = sqrt(dot(l, l, 1)); % length of the dipole
[ lx, ly, lz ] = uncat(1, l./ll); % dipole direction vector
[ ex ey ez ] = dfield(freq,eps0,mu0,dx,dy,dz,lx,ly,lz,ll);
dzi = -rsrc(3)-robs(3)+2*lay.z2;  % image position
[ lxi, lyi ] = deal(-lx, -ly); % flipped in transverse direction
[ exi eyi ezi ] = dfield(freq,eps0,mu0,dx,dy,dzi,lxi,lyi,lz,ll);
ef2 = [ ex+exi; ey+eyi; ez+ezi ];

assertEquals(ef2, ef1, 1e-10);

% Second case - perfectly conducting plane on top
lay.z3   = 3e-3;
lay.z2   = 1e-3;
lay.eps3 = eps0-j*Inf;
lay.eps2 = eps0;
lay.eps1 = eps0;

gi = mkimages(lay);

ef1 = dfield_dcim(lay,gi,robs,rsrc,l);

[ dx, dy, dz ] = uncat(1, rsrc-robs);
ll = sqrt(dot(l, l, 1)); % length of the dipole
[ lx, ly, lz ] = uncat(1, l./ll); % dipole direction vector
[ ex ey ez ] = dfield(freq,eps0,mu0,dx,dy,dz,lx,ly,lz,ll);
dzi = -rsrc(3)-robs(3)+2*lay.z3;  % image position
[ lxi, lyi ] = deal(-lx, -ly); % flipped in transverse direction
[ exi eyi ezi ] = dfield(freq,eps0,mu0,dx,dy,dzi,lxi,lyi,lz,ll);
ef2 = [ ex+exi; ey+eyi; ez+ezi ];

assertEquals(ef2, ef1, 1e-10);

% Consider adding the case of planes both on top and bottom. In this case
% there are two series of images - the source dipole is mirrored by the
% first plane, then this image is mirrored by the second plane, then by
% the first one again etc.

