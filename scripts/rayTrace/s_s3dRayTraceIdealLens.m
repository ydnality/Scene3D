%% Ray-tracing for an ideal lens
%
% An ideal thin lens is defined as one that obey's the thin lens equation
% (1/s1 + 1/s2 = 1/f).  Therefore, Snell's Law is not used in this script
% at all.  Instead, the direction that the rays bend are determined by the
% thin lens equation.  At each point source, a "center ray" is shot at the
% center of the lens.  The intersection of this ray and the focal-plane as
% defined by the thin lens equation determines the point of focus of the
% lens.  All other rays that are shot at the edge of the aperture will then
% intersect this ideal focal point.  
% All units are in mm.
%
%  TODO:
%   We can't use sensor here and sensorCreate in same way.
%   Somehow we should integrate this sensor stuff with ISET.
%   I don't see the lens part here.  Need more comments about what this is
%   doing, and what it tests.
%
% See also: s_3dRayTrace*.m
%
% AL Vistalab 2014

%
ieInit

%% point sources
% declare point sources in world space.  The camera is usually at [0 0 0],
% and pointing towards -z.  We are using a right-handed coordinate system.

% [XGrid YGrid] = meshgrid(-4000:1000:4000,-4000:1000:4000);
[XGrid YGrid] = meshgrid(-2000:1000:2000,-2000:1000:2000);  %large
% distance points
% [XGrid YGrid] = meshgrid(-10:100:10,-10:10:10);  %small distance points
pointSources = [XGrid(:) YGrid(:) ones(size(XGrid(:))) * -20000];
% pointSources = [XGrid(:) YGrid(:) ones(size(XGrid(:))) * -100];
% pointSources = [0 0 -20000];
wave = 400:10:700;
rtType = 'ideal';  %ray tracing using the ideal lens model  
% rtType = 'realistic';

%% camera and film properties 

% Build a sensor (film) object
% Position, size,  wave, waveConversion, resolution
film = filmC('position', [0 0 50],'size', [7 7], 'wave', wave);   

diffractionEnabled = true;
% lensObject('ideal')
% apertureRadiusMM = .1; %diffraction case
apertureMiddleD = .1;
sensorDistance = 50;  %mm
apertureSamples = 100;

% lens = lensIdealObject(apertureRadiusMM, sensorDistance, [0 0 0], diffractionEnabled);
lens = lensC('apertureMiddleD', apertureMiddleD,  ...
    'diffractionEnabled', diffractionEnabled, 'wave', wave, 'apertureSample', [apertureSamples apertureSamples]);


%% loop through all point sources

for curInd = 1:size(pointSources, 1);
    
    %calculate the origin and direction of the rays
    rays = lens.rtSourceToEntrance(pointSources(curInd, :), [], [], rtType);
    
    %duplicate the existing rays, and creates one for each
    %wavelength
    rays.expandWavelengths(film.wave);

    %lens intersection and raytrace
    lens.rtThroughLens(rays, [], rtType);

    % intersect with "film" and add to film
    rays.recordOnFilm(film);
end

%% Assign to optical image
oi = oiCreate;
oi = initDefaultSpectrum(oi);
oi = oiSet(oi,'photons',film.image);

% Set the optics parameters too
%
optics = oiGet(oi,'optics');
optics = opticsSet(optics,'focal length',lens.focalLength/1000);
optics = opticsSet(optics,'fnumber', lens.focalLength/(2*apertureMiddleD));
oi = oiSet(oi,'optics',optics);

% Opposite over adjacent is the tan of half the angle ...
% Everything is mm
% hfov = rad2deg(2*atan2(apertureRadiusMM,lens.focalLength));
hfov = rad2deg(2*atan2(film.size(1)/2,lens.focalLength));
oi = oiSet(oi,'hfov', hfov);

vcAddObject(oi); oiWindow;

%%
