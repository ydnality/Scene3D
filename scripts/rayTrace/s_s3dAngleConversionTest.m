%% Testing of the spherical angle conversion function in the rayC class
%
% This script renders a PSF and in the process tests the spherical angle
% conversion function in the rayC class
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
% film = pbrtFilmC([0 0 51.2821	],[.2/sqrt(2) .2/sqrt(2)], 400:10:700, [(400:10:700)' (1:31)'], [50 50 31]);

% Declare film
filmPosition = [0 0 51.2821	];
filmSize = [.2/sqrt(2) .2/sqrt(2)];
wave = 400:10:700;
resolution =  [50 50 length(wave)];
film = pbrtFilmC('position', filmPosition, 'size', filmSize, 'wave', wave, 'resolution', resolution);

% Declare Lens
diffractionEnabled = false;   %diffraction causes imaginary directions!! TODO:  INVESTIGATE THIS!!
apertureDiameterMM = 2.2727;  %f/22 
% apertureDiameterMM = 3.1250;  %f/16
% apertureDiameterMM = 4.5455;  %f/11

fLength = 50;
apertureSamples = [51 51];
name = 'idealLensTest';
type = 'idealLens';
jitterFlag = true;
debugLines = 100;
lens = lensC('name', name, 'type', type, 'focalLength', fLength, 'diffractionEnabled', diffractionEnabled, 'wave', wave, 'aperturesample', apertureSamples);

%% Loop through all point sources and Render PSF
vcNewGraphWin; 

curInd = 1;

%% calculate the origin and direction of the rays
disp('-----trace source to lens-----');
tic
rtType = 'ideal';
rays = lens. rtSourceToEntrance(pointSources(curInd, :), false, jitterFlag, rtType)
toc

%%  look at angles and analyze them
sphereAngles = rays.get('sphericalAngles');

figure; hist(sphereAngles(:,1));  %should be uniform
%should have more at edges (more data points at perimeter) ... approximately linear
figure; hist(sphereAngles(:,2));  

%% use the projection form of nagles

projAngles = rays.get('projectedAngles');
hist(projAngles(:,1)); 
hist(projAngles(:,2)); 


%% plot phase space

rays.plotPhaseSpace();

%% duplicate the existing rays, and creates one for each wavelength
disp('-----expand wavelenghts-----');
tic
rays.expandWavelengths(film.wave);
toc
disp('-----rays trace through lens-----');
tic

%lens intersection and raytrace
lens.rtThroughLens(rays, debugLines, rtType);
toc

%% plot phase space
rays.plotPhaseSpace();


%%  intersect with "film" and add to film
disp('-----record on film-----');
tic
rays.recordOnFilm(film);
toc


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
