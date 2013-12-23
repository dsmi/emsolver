function intg = ilay_c(lay, gi, r, robs, qN)
% intg = ilay_c(lay, gi, r, robs, qN)
%
% Evaluates integral of the correction term of the Greens function in layered
% media.
%

imga = { gi.Cf1.a gi.Cf2.a gi.Cf3.a gi.Cf4.a };
imgb = { gi.Cf1.b gi.Cf2.b gi.Cf3.b gi.Cf4.b };

intg = integ_lay(@integ_fp, lay, 0, imga, imgb, r, robs, qN);
