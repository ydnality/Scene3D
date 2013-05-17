%% Scratch code for checking the math of known intensity ratio to distances
%
%
%

% Assume flash 1 is a point at the origin
% Assume flash 2 is a point at (0,0,-f)
% Suppose d1 is the distance to flash 1 and d2 is the distance to flash 2.

% Consider a point, p, at 
p = [3,4,5];

% We can express the location of point p in terms of its distance from
% flash 1 and flash 2, and also in terms of the azimuth and elevation of
% the point with respect to the origin.  (We should make a picture on
% google docs to illustrate the situation).

% For any point in the sensor plane, (x,y,-1)
p = [1,1,1]
d1 = norm(p);

% Suppose that camera plane is at (0,0,1) and the image model is a pinhole
% camera where the pinhole is also at the origin, like the flash.
% Then we should be able to derive the angles, alpha and phi for every
% pixel in the camera plane.  

alpha = asin(p(2)/d1); 
phi   = asin(p(3)/(d1*cos(alpha)));
%
% This is a check
% q(1) = d1*cos(alpha)*cos(phi);
% q(2) = d1*sin(alpha);    
% q(3) = d1*cos(alpha)*sin(phi);
% q

% If someone gives us an intensity ratio at a point for the two images,
% call it, iRatio, for flash 1 over flash 2, we believe that is equal to
% (d2/d1)^2.   We can estimate the depth map, which is d1, as

f =  1;   % f is one unit behind the origin
iRatio = (1:.01:5);
d1 = zeros(size(iRatio));
% iRatio(1) = 1;
for ii=1:length(iRatio)
    radical = sqrt(4*cos(alpha).^2*sin(phi).^2.*f.^2  - 4*f^2.*(1 - iRatio(ii)));
    d1(ii) = (2*f^2)/(-2*cos(alpha).*sin(phi).*f + radical);
end

% We have a problem when phi is negative.  Otherwise, looking pretty good.

vcNewGraphWin;
loglog(iRatio,d1)
xlabel('iRatio')
ylabel('Estimated distance')
grid on




%% testing the equation using known values

% assume flash 1 at the origin
% flash 2 is at (0,0,-1)

% consider a point at (0,0,1)
% here, at this point, alpha = 0, phi= pi/2

% according to r^2 fall off, picture 1 should have 4 times the radiance as
% picture 2.  1/(r^2 = 1) vs. 1/(r^2 = 4)  
% thus, iRatio for this point should be 4, and d1 = 1


testIRatio = 4;
f = 1;
phi = pi/2;
alpha = 0;
radical = sqrt(4*cos(alpha).^2*sin(phi).^2.*f.^2 - 4*f^2.*(1 - testIRatio));
d1Test = (2*f^2)/(-2*cos(alpha).*sin(phi).*f + radical)

% this checks out, d1Test = 1, as expected!


% flash 1 at origin -> intensity ~ 1/2^2
% flash 2 at (0, 0, -4) -> intensity ~ 1/6^2
% point at (0, 0, 2)
% alpha = 0, phi = pi/2

% picture 1 should have 9 times the radiance of picture 2
% thus, iRatio for this point should be 9, and d1 = 2

testIRatio = 9;
f = 4;
phi = pi/2;
alpha = 0;
radical = sqrt(4*cos(alpha).^2*sin(phi).^2.*f.^2 - 4*f^2.*(1 - testIRatio));
d1Test = (2*f^2)/(-2*cos(alpha).*sin(phi).*f + radical)

% this checks out, d1Test = 2, as expected!



% flash 1 at origin -> intensity ~ 1/50
% flash 2 at (0, 0, -6) -> intensity ~ 1/ 146
% point at (3, 4, 5)
% alpha = 0, phi = pi/2

% picture 1 should have 2.92 times the radiance of picture 2
% thus, iRatio for this point should be 2.92, and d1 = 7.0711

p = [3 4 5]; d1 = 7.0711;
testIRatio = 2.92;
f = 6;
alpha = asin(p(2)/d1)
phi   = asin(p(3)/(d1*cos(alpha)))

radical = sqrt(4*cos(alpha).^2*sin(phi).^2.*f.^2 - 4*f^2.*(1 - testIRatio));
d1Test = (2*f^2)/(-2*cos(alpha).*sin(phi).*f + radical)

% this checks out, d1Test = 7.0711, as expected!

% flash 1 at origin -> intensity ~ 1/50
% flash 2 at (0, 0, -6) -> intensity ~ 1/ 146
% point at (3, 4, 5)
% alpha = 0, phi = pi/2

% picture 1 should have 2.92 times the radiance of picture 2
% thus, iRatio for this point should be 2.92, and d1 = 7.0711

p = [-3 -4 5]; d1 = 7.0711;
testIRatio = 2.92;
f = 6;
alpha = asin(p(2)/d1)
phi   = asin(p(3)/(d1*cos(alpha)))

radical = sqrt(4*cos(alpha).^2*sin(phi).^2.*f.^2 - 4*f^2.*(1 - testIRatio));
d1Test = (2*f^2)/(-2*cos(alpha).*sin(phi).*f + radical)

% this checks out, d1Test = 7.0711, as expected!