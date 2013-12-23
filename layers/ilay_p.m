function intg = ilay_p(lay, gi, r, robs, qN)
% intg = ilay_p(lay, gi, r, robs, qN)
%
% Evaluates integral of the Greens function for scalar potential in layered
% media. The semantics of this function is similar to integ_p which evaluates
% integral of the free-space Green's function, this function is supposed
% to be called instead of integ_p when modeling the layered dielectric.
% The spatial domain Green's function for layered media is approximated with
% complex images (see mkimages), so under the hood this function calls
% integ_p for each of the images.
%

m = 2*(-j*lay.freq*eps0);

imga = { m*gi.Kf1.a m*gi.Kf2.a m*gi.Kf3.a m*gi.Kf4.a };
imgb = {   gi.Kf1.b   gi.Kf2.b   gi.Kf3.b   gi.Kf4.b };

a0 = eps0/lay.eps2;
intg = integ_lay(@integ_p, lay, a0, imga, imgb, r, robs, qN);
