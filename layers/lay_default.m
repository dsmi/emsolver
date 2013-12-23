function lay = lay_default
% lay = lay_default
%
% Returns the structure which describes the layered media parameters
% initialized by default.
% The current implementation of the Green's function calculator supports
% up to three layers, conductors are buried into the middle layer.
% The layers are numbered from top to bottom along the Z axis positive
% direction.
%
% The strucure contains the following fields:
%   freq     Angular frequency
%   z2       Z-coordinate of the interface between the first and
%            the second layers.
%   z3       Z-coordinate of the interface between the second and
%            the third layers. Layers are numbered from top to bottom
%            along the Z axis positive direction, so lay_z3 > lay_z2.
%   eps1 eps2 eps3  Electric permittivity of the layers.
%

lay.freq = 1e10;
lay.z2 = 0.01;
lay.z3 = 0.02;
lay.eps1 = eps0*5;
lay.eps2 = eps0*2;
lay.eps3 = eps0*12;
