% 4th Order Wavefront Aberration


function [W]=pa4thWaveAber(k,l,ro,theta)

%INPUT
%k: order of radius
%l: order of theta
%ro: radius
%theta: theta

%OUTPUT
%W: Term of 4th order Wavefront aberration

W=ro.^k.*cos(theta).^l;