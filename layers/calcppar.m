function [ T01, T02 ] = calcppar(lay)
% [ T01, T02 ] = calcppar(lay)
%
% This function calculates the integration path paramaters T01 and T02
% based on the layers stack declared in the input structure (see lay_default.m).
% For the details on the integration path used and meaning of the T01 and T02
% parmeters please see mkipath1 and mkipath2
%

lay_k = calc_lay_k(lay);
iobs = 2; % Observation layer hardcoded
ki = lay_k(iobs);

% Remove infinite entries - to handle ideal grounds properly
lay_k_noinf = lay_k(find(abs(lay_k) < 1e50));

max_k_ratio = max(real(lay_k_noinf/ki));
T02 = sqrt((1.5*max_k_ratio).^2-1);
T01 = sqrt((100*max_k_ratio).^2-1)-T02;
