function intg = ilay_fp(lay, gi, r, robs, qN)
% intg = ilay_fp(lay, gi, r, robs, qN)
%
% Evaluates integral of the Greens function for vector potential in layered
% media. The semantics of this function is similar to integ_fp which evaluates
% integral of the free-space Green's function, this function is supposed
% to be called instead of integ_fp when modeling the layered dielectric.
% The spatial domain Green's function for layered media is approximated with
% complex images (see mkimages), so under the hood this function calls
% integ_fp for each of the images.
%

m = 2/(j*lay.freq*mu0);

imga = { m*gi.Gvv1.a m*gi.Gvv2.a m*gi.Gvv3.a m*gi.Gvv4.a };
imgb = {   gi.Gvv1.b   gi.Gvv2.b   gi.Gvv3.b   gi.Gvv4.b };

igvv = integ_lay(@integ_fp, lay, 1, imga, imgb, r, robs, qN);
igvv(:,:,3) = 0;


imga = { m*gi.Gzz1.a m*gi.Gzz2.a m*gi.Gzz3.a m*gi.Gzz4.a };
imgb = {   gi.Gzz1.b   gi.Gzz2.b   gi.Gzz3.b   gi.Gzz4.b };

igzz = integ_lay(@integ_fp, lay, 1, imga, imgb, r, robs, qN);
igzz(:,:,1:2) = 0;

intg = igvv+igzz;

% Notice the -1/j multiplier - differentiation by dx/dy in spatial
% domain corresponds to multiplication by -j*kx / -j*ky!
m = (-1/j)*2/(j*lay.freq*mu0);

imga = { m*gi.Gzu1.a m*gi.Gzu2.a m*gi.Gzu3.a m*gi.Gzu4.a };
imgb = {   gi.Gzu1.b   gi.Gzu2.b   gi.Gzu3.b   gi.Gzu4.b };

% Calculate div(fxy)
iaux = integ_aux(r, robs, qN);
nxy = iaux.n;
nxy(:,:,3) = 0;
divfxy = 2-sum(nxy.*nxy, 3);

% div(fxy)*\int g dS
div_intg = divfxy.*integ_lay(@integ_p, lay, 0, imga, imgb, r, robs, qN);
% \int div(fxy * g) dS
intg_div = integ_lay(@integ_fxy, lay, 0, imga, imgb, r, robs, qN);

intg(:,:,3) = intg(:,:,3) + intg_div - repmat(div_intg, 1, 3);
