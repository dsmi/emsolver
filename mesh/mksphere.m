function [tri, x, y, z] = mksphere(n)
% [tri, x, y, z] = mksphere(n)
%
% Makes a 3D sphere. Wrapper of 3-rd party BuildSphere.m

[p,t]=BuildSphere(n);
tri = [ t(:,2) t(:,1) t(:,3) ]; % To get the correct orientation.
x = p(:,1).';
y = p(:,2).';
z = p(:,3).';
