% Tests the calculation of spherical lens properties
% Calculates the radius of curvature needed for a 2 element lens lens with focal length f,
% index of refraction n, and lens element selaration d.  The 2 lens
% elements will both have a radius of curvature R (one will be R, the other
% will be -R) 

n = 1.67;   %Index of refraction (mm)
f = .259 %.28; %.259; %.23739;     %Focal length (mm)
d = .01;    %Lens element separation (mm)

f = f + d/2

oneOverR = (-2 * (n-1) + sqrt(4 * (n-1)^2 + 4 * (n-1)^2 * d/(n * f)))/(2 * (n-1)^2 * d/n)
R = 1/oneOverR

check = (n-1) * (2/R + (n-1) * d/(n*R^2))


%%  we need to check that the lenses are not too curved for the lenslet size

xInt = -(2 * R - d)/2  %this is the x intercept point of where the 2 lens elements intersect
yInt = sqrt(R^2 -xInt^2)

maxHeight = yInt * 2    % in mm