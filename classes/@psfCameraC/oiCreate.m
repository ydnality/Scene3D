function oi = oiCreate(obj)
% Create an ISET optical image from the psfCamera object
%
% The PSF should be calculated prior to this if you would like to see the
% point spread in the oi, which is the typical usage
%
% AL/BW Vistasoft Team, Copyright 2014

% Create an optical image from the camera (film) image data.
oi = oiCreate;
oi = initDefaultSpectrum(oi);
oi = oiSet(oi,'wave', obj.film.wave);
oi = oiSet(oi,'photons',obj.film.image);

% The photon numbers do not yet have meaning.  This is a hack,
% that should get removed some day, to give the photon numbers
% some reasonable level.
oi = oiAdjustIlluminance(oi,1);

% Set focal length in meters
oi = oiSet(oi,'optics focal length',obj.lens.focalLength/1000);

% This isn't exactly the fnumber.  Do we have the aperture in
% there?  Ask MP what to use for the multicomponent system.
% This is just a hack to get something in there
% Maybe this should be obj.lens.apertureMiddleD?
fN = obj.lens.focalLength/obj.lens.surfaceArray(1).apertureD;
oi = oiSet(oi,'optics fnumber',fN);

% Estimate the horizontal field of view
hfov = rad2deg(2*atan2(obj.film.size(1)/2,obj.lens.focalLength));
oi = oiSet(oi,'hfov', hfov);

% Set the name based on the distance of the sensor from the
% final surface.  But maybe the obj has a name, and we should
% use that?  Or the film has a name?
temp = obj.film.position;
filmDistance = temp(3);
oi = oiSet(oi, 'name', ['filmDistance: ' num2str(filmDistance)]);

end