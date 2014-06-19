function [xp,yp] = circlePoints(center,r, aStep)
% Draw a circle of radius r
%
% center coordinates of the center of the circle: default (0,0)
% r is the radius of the circle:             default 1
% aStep is the angular step in radians:        default is 0.01
%
% Example:
%    [x,y] = circlePoints; plot(x,y,'.'); axis equal
%    [x,y] = circlePoints([2,2],3,0.1); ; plot(x,y,'.'); axis equal
%
% AL/BW Vistasoft team 2014

if ieNotDefined('center'), center = [0,0]; end
if ieNotDefined('r'), r = 1; end
if ieNotDefined('ang'), aStep = 0.01; end

% Angles and points
ang = 0:aStep:2*pi; 
xp  = r*cos(ang);
yp  = r*sin(ang);

% Shift center
xp = xp + center(1); yp = yp + center(2);

end