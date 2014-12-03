function pts = ptsS2P(pts)
% Convert 3D point (space) into a polar represetation
%
%   ptsPolar = ptsS2P(pts)
%
% AL Vistasoft 2014

% For the moments, pts is just a 3-vector
currentDepth = sqrt(pts(1)^2 + pts(2)^2 + pts(3)^2);

% We are worried about the negative sign here ... and whether z is always
% negative or not ...
currentAngle = atan(pts(2)/(-pts(3))); 

pts = [0 currentAngle currentDepth];

end