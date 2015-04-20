function oi = s3dFixOi(oi, inParams)
% Adjust the oi for consistency with pbrt calculations
%
%   oi = s3dFixOi(oi, [pbrtO or focalLength])
% 
% This should not be needed some day.  It is used in s3dRenderOI for now.
%
% pbrtO:  In the preferred use, we send in a pbrt object and derive the
% various parameters for the oi
%   It is also possible to just send in the focalLength: the theoretical
%   focal length of the lens.  This becomes the sensor distance.  But this
%   is an unpreferred hack that assumes the lens is focused at infinity.
%
% See also: pbrt2oi
%
% Vistalab 2015

error('Deprecated.  See pbrt2oi');

end

% Scale the values ... not sure why this level.
oi = oiSet(oi, 'photons', oiGet(oi,'photons') * 10^13);  %some normalization issues

% We need to figure out the field of view
if isequal(class(inParams),'pbrtObject')
    % Compute the horizontal field of view
    fdiag = inParams.camera.lens.filmDiag;
    dist = inParams.camera.lens.filmDistance;
    x = inParams.camera.film.xresolution;
    y = inParams.camera.film.yresolution;
    d = sqrt(x^2+y^2);
    fwidth = (fdiag/d)*x;
    fov = 2*atan2d(fwidth/2,dist);
elseif isnumeric(inParams)
    % This should never happen!  Alert the user that we are completely
    % making this up.
    dist = inParams;  % Treat the focal length as the location of the sensor
    fwidth = 30;      %  Thirty five millimeter diagonal, 
    fov = 2*atan2d(fwidth/2,dist);
    fprintf('Arbitrarily setting the film width to %.2f mm',fwidth)
else
    error('Unknown input parameter type');
end

oi = oiSet(oi,'fov', fov);

    
end