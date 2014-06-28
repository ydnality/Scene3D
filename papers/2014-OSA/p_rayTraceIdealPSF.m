%% 2014 OSA Conference: Ray-tracing for an ideal lens
%
% This script renders the PSF for an ideal thin lens under various point
% source positions.  The PSF is displayed as an optical image (no sensor yet).
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
%
% See also: s_3dRayTrace*.m
%
% For Figure PSF comparison
%
% AL Vistalab 2014

%% Initialize Iset
s_initISET

%% Declare point sources
% tic
% declare point sources in world space.  The camera is usually at [0 0 0],
% and pointing towards -z.  We are using a right-handed coordinate system.

% [XGrid YGrid] = meshgrid(-4000:1000:4000,-4000:1000:4000);
[XGrid YGrid] = meshgrid(-2000:1000:2000,-2000:1000:2000);
pointSources = [0 0 -2000];
% pointSources = [0 0 -20000];


%% Declare camera and film properties

% Build a sensor (film) object
% Position, size,  wave, waveConversion, resolution
% film = pbrtFilmObject([0 0 51.2821	],[.2/sqrt(2) .2/sqrt(2)], 400:10:700, [(400:10:700)' (1:31)'], [50 50 31]);

% Declare film
filmPosition = [0 0 51.2821	];
filmSize = [.2/sqrt(2) .2/sqrt(2)];
wave = 400:10:700;
resolution =  [50 50 length(wave)];
film = pbrtFilmObject('position', filmPosition, 'size', filmSize, 'wave', wave, 'resolution', resolution);

% Declare Lens
diffractionEnabled = true;
apertureDiameterMM = 2.2727;  %f/22 
% apertureDiameterMM = 3.1250;  %f/16
% apertureDiameterMM = 4.5455;  %f/11

fLength = 50;
apertureSamples = [51 51];
name = 'idealLensTest';
type = 'idealLens';
jitterFlag = true;
debugLines = 100;
lens = lensMEObject('name', name, 'type', type, 'focalLength', fLength, 'diffractionEnabled', diffractionEnabled, 'wave', wave, 'aperturesample', apertureSamples);

%% Loop through all point sources and Render PSF
vcNewGraphWin; 

for curInd = 1:size(pointSources, 1);

    %calculate the origin and direction of the rays
    disp('-----trace source to lens-----');
    tic
    rtType = 'ideal';
    rays = lens. rtSourceToEntrance(pointSources(curInd, :), false, jitterFlag, rtType)
    toc
    
    %duplicate the existing rays, and creates one for each
    %wavelength
    disp('-----expand wavelenghts-----');
    tic
    rays.expandWavelengths(film.wave);
    toc
    disp('-----rays trace through lens-----');
    tic
    
    %lens intersection and raytrace
    lens.rtThroughLens(rays, debugLines, rtType);
    toc
    
    % intersect with "film" and add to film
    disp('-----record on film-----');
    tic
    rays.recordOnFilm(film);
    toc
end

%% Assign to optical image

oi = oiCreate;
oi = initDefaultSpectrum(oi);
oi = oiSet(oi,'photons',film.image);

% Set the optics parameters too
optics = oiGet(oi,'optics');
optics = opticsSet(optics,'focal length',lens.focalLength/1000);
optics = opticsSet(optics,'fnumber', lens.focalLength/(apertureDiameterMM));
oi = oiSet(oi,'optics',optics);

% Opposite over adjacent is the tan of half the angle ...
% Everything is mm
% hfov = rad2deg(2*atan2(apertureRadiusMM,lens.focalLength));
hfov = rad2deg(2*atan2(film.size(1)/2,lens.focalLength));
oi = oiSet(oi,'hfov', hfov);

vcAddObject(oi); oiWindow;
