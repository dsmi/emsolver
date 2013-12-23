function [ kr, kz, t ] = mkipath1(lay, T01, T02, ns)
% function [ kr, kz, t ] = mkipath1(lay, T01, T02, ns)
%
% Layered Green's funtions calculator employs the two-level approximation
% scheme - the integration/approximation path in kr plane consists of
% two parts; first one is along real positive kr which starts from
% krmax2 and goes to krmax1, the second part is curved and goes from
% 0 to krmax2.
% This function creates the first part of the integration path.
%   Inputs:
%     lay - structure describing the layered media parameters.
%     T01 - Integration path parameter, is calculated based on
%           krmax1.
%     T02 - Integration path parameter, is calculated based on
%           krmax2.
%     ns  - number of the samples in the path
%   Outputs:
%     kr  - radial component of the wavevector, column vector of
%           length ns.
%     kz  - vertical (perpendicular to the stratification) component of
%			the wavevector, ns-by-num_of_layers array.
%     t   - real variable which is uniformly sampled over the range
%           [ 0, T01 ]. The path is defined by mapping from t variable to
%           the complex kr plane. This variable is neccessary because
%           GPOF requires uniform sampling over a real variable. Column
%           vector of length ns.
%

lay_k = calc_lay_k(lay);
iobs = 2;
ki = lay_k(iobs);

t=linspace(0, T01, ns).';
kr2 = 1 + (T02+t).^2;
kr = ki*sqrt(kr2);
kz = calc_kz(lay, kr);
