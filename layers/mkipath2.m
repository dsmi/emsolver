function [ kr, kz, t ] = mkipath2(lay, T02, ns)
% function [ kr, kz, t ] = mkipath2(lay, T02, ns)
%
% Layered Green's funtions calculator employs the two-level approximation
% scheme - the integration/approximation path in kr plane consists of
% two parts; first one is along real positive kr which starts from
% krmax2 and goes to krmax1, the second part is curved and goes from
% 0 to krmax2.
% This function creates the second part of the integration path.
%   Inputs:
%     lay - structure describing the layered media parameters.
%     T02 - Integration path parameter, is calculated based on
%           krmax2.
%     ns  - number of the samples in the path
%   Outputs:
%     kr  - radial component of the wavevector, column vector of
%           length ns.
%     kz  - vertical (perpendicular to the stratification) component of
%			the wavevector, ns-by-num_of_layers array.
%     t   - real variable which is uniformly sampled over the range
%           [ 0, T02 ]. The path is defined by mapping from t variable to
%           the complex kr plane. This variable is necessary because
%           GPOF requires uniform sampling over a real variable. Column
%           vector of length ns.
%

lay_k = calc_lay_k(lay);
iobs = 2;
ki = lay_k(iobs);

t=linspace(0, T02, ns).';
kr2 = 1-(-j*t+(1-t/T02)).^2;
kr = ki*sqrt(kr2);
lkr = repmat(lay_k, ns, 1);
kz = calc_kz(lay, kr);

