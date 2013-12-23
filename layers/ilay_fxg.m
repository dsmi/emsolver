function intg = ilay_fxg(lay, gi, r, robs, qN)
% intg = ilay_fxg(lay, gi, r, robs, qN)
%
% Evaluates integral of the curl of Greens function for electric vector
% potential in layered media. The semantics of this function is similar
% to integ_fxg which evaluates integral of the curl of free-space Green's
% function, this function is supposed to be called instead of integ_fxg
% when modeling the layered dielectric.
% The spatial domain Green's function for layered media is approximated with
% complex images (see mkimages), so under the hood this function calls
% integ_fxg for each of the images.
%

m = -j*2/(lay.eps2*lay.freq);

imga = { m*gi.Fvv1.a m*gi.Fvv2.a m*gi.Fvv3.a m*gi.Fvv4.a };
imgb = {   gi.Fvv1.b   gi.Fvv2.b   gi.Fvv3.b   gi.Fvv4.b };

Mvv = diag([ 1 1 0 ]);
integ_fxg_vv = @(k, r, robs, qN) integ_tfxg(k, r, robs, Mvv, qN);
ifvv = integ_lay(integ_fxg_vv, lay, 1, imga, imgb, r, robs, qN);

imga = { m*gi.Fzz1.a m*gi.Fzz2.a m*gi.Fzz3.a m*gi.Fzz4.a };
imgb = {   gi.Fzz1.b   gi.Fzz2.b   gi.Fzz3.b   gi.Fzz4.b };

Mzz = diag([ 0 0 1 ]);
integ_fxg_zz = @(k, r, robs, qN) integ_tfxg(k, r, robs, Mzz, qN);

ifzz = integ_lay(integ_fxg_zz, lay, 1, imga, imgb, r, robs, qN);

intg = ifvv+ifzz;

% Notice the -1/j multiplier - differentiation by dx/dy in spatial
% domain corresponds to multiplication by -j*kx / -j*ky!
m = (-1/j)*m;

imga = { m*gi.Fzu1.a m*gi.Fzu2.a m*gi.Fzu3.a m*gi.Fzu4.a };
imgb = {   gi.Fzu1.b   gi.Fzu2.b   gi.Fzu3.b   gi.Fzu4.b };

intg_fzu = integ_lay(@integ_cfxy, lay, 0, imga, imgb, r, robs, qN);

% Minus, because curl(f*g) = -cross(f, grad(g))
intg = intg - intg_fzu;
