function [ ex, ey, ez, hx, hy, hz ] = dfield(freq,eps,mu,x,y,z,dx,dy,dz,dl)
% [ ex, ey, ez, hx, hy, hz ] = dfield(freq,eps,mu,x,y,z,dx,dy,dz,dl)
%
% Compute the fields of the electric current dipole.
% freq         - angular frequency
% eps, mu      - parameters of the media
% R(x,y,z)     - centers of the dipole(s)
% R(dx,dy,dz)  - orientation of the dipoles
% dl           - length of the dipoles
%
% Observation point is located at (0,0,0).

k = freq * sqrt(eps * mu);

% Scalar potential of the charge at the beginning of the dipole.
% Q1 = -I/j*freq (integral of delta*div(I) = I)
dp1 = -1/(j*freq)*diffp(k, -(x-dx.*dl/2), -(y-dy.*dl/2), -(z-dz.*dl/2));
[ dpx1 dpy1 dpz1 ] = uncat(ndims(dp1), shiftdim(dp1, 1));

% Scalar potential of the charge at the end of the dipole.
% Q2 = I/j*freq (integral of delta*div(I) = -I)
dp2 =  1/(j*freq)*diffp(k, -(x+dx.*dl/2), -(y+dy.*dl/2), -(z+dz.*dl/2));
[ dpx2 dpy2 dpz2 ] = uncat(ndims(dp2), shiftdim(dp2, 1));

% Vector potential of the current
A = mu/(4*pi)*calcp(k, x, y, z).*dl; % magnitude
Ax = A.*dx;
Ay = A.*dy;
Az = A.*dz;

% Compute the electric field from the vector and scalar potentials
ex = -j*freq*Ax - ( dpx1 + dpx2 )/(eps*4*pi);
ey = -j*freq*Ay - ( dpy1 + dpy2 )/(eps*4*pi);
ez = -j*freq*Az - ( dpz1 + dpz2 )/(eps*4*pi);

% Next, compute the magnetic field as H=(1/mu)*curl(A)
% curlA = i*(dAz/dy-dAy/dz)+j*(dAx/dz-dAz/dx)+k*(dAy/dx-dAx/dy)
% dAx/dx = dA/dx*x, dAy/dx = dA/dx*y, dAz/dx = dA/dx*z, etc
% AxD = i*(Ay*Dz-Az*Dy)+j*(Az*Dx-Ax*Dz)+k*(Ax*Dy-Ay*Dz)
dA = mu/(4*pi)*diffp(k, x, y, z).*repmat(shiftdim(dl, -1), 3, 1);
d = cat(1, shiftdim(dx, -1), shiftdim(dy, -1), shiftdim(dz, -1));
curlA = cross(dA, d, 1);
H = shiftdim(curlA,1)/mu;
[ hx hy hz ] = uncat(ndims(H), H);
