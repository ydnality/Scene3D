function d = rayDirection(oPoints,ePoints)
% Compute direction vector between origin and end points of a ray
%
% AL/BW (c) Vistasoft Team, 2014

d = ePoints - oPoints;
d = normvec(d,'dim',2);

end