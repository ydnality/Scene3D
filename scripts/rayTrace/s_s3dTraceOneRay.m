%% Testing of the spherical angle conversion function in the rayC class
%
% This script ray-traces a single ray from the point source, through the
% lens, to the sensor.  
%
% The resulting rays object will contain that one ray, but at the exit
% pupil.  
%
% AL Vistalab 2014

%% Initialize Iset
s_initISET

%% Declare single ray direction from point source

wave = 400:10:700;  %wavelength samples that all components can accomodate
rayDirection = [0 .1 1]; %direction for the single ray we want to trace
rayDirection = normvec(rayDirection);
rayWavelength = 500;

waveIndex = find(wave == rayWavelength)


%% Declare ray-trace type

rtType = 'realistic';  %ideal/realistic
debugLines = 50;
%% Declare point sources
% declare point sources in world space.  The camera is usually at [0 0 0],
% and pointing towards -z.  We are using a right-handed coordinate system.

% Use the psCreate thing ...
%pointSources = [0 10 -2000];
pointSources = [0 0 -50];

%% Declare camera properties

% Declare film
% filmPosition = [0 0 51.2821	];  % Good for 2Elens
%filmPosition = [0 0 37.4];  % Good for dgauss.50mm.  True focal about 37.3mm
filmPosition = [0 0 80];  % Good for dgauss.50mm.  True focal about 37.3mm

% Build a sensor (film) object
% filmSize = [.2/sqrt(2) .2/sqrt(2)];
%filmDiag = 1;  % Millimeters
filmDiag = 40;  % Millimeters
filmSize = [filmDiag/sqrt(2) filmDiag/sqrt(2)];

resolution =  [300 300 length(wave)];
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

thinlens = lensC('name', name, 'type', type, 'focalLength', fLength, 'diffractionEnabled', diffractionEnabled, 'wave', wave, 'aperturesample', apertureSamples);

lensFile = fullfile(s3dRootPath, 'data', 'lens', '2ElLens.mat');
import = load(lensFile,'lens');
thickLens = import.lens;
thickLens.apertureMiddleD = 10;

lensFile = fullfile(s3dRootPath, 'data', 'lens', 'dgauss.50mm.mat');
import = load(lensFile,'lens');
multiLens = import.lens;

% lens = thickLens;
lens = multiLens;
lens.set('wave', wave);  %TODO: right now only 400:100:700 works.  Fix this.

%% calculate the origin and direction of the rays
curInd = 1;
disp('-----trace source to lens-----');
tic

%rays = lens.rtSourceToEntrance(pointSources(curInd, :), false, jitterFlag, rtType)
rays = rayC('origin', pointSources(curInd, :), 'direction', rayDirection, ...
    'wave', wave, 'waveIndex', 1);
toc

%% use the projection form of angles and plot phase space
projAngles = rays.get('projectedAngles');

% plot phase space
rays.plotPhaseSpace();

%% duplicate the existing rays, and creates one for each wavelength
disp('-----expand wavelenghts-----');
% tic
% rays.expandWavelengths(film.wave);
% toc
disp('-----rays trace through lens-----');
tic

%lens intersection and raytrace
lens.rtThroughLens(rays, debugLines, rtType);
toc

% plot phase space
rays.plotPhaseSpace();

%%  intersect with "film" and add to film

disp('-----record on film-----');
tic
rays.recordOnFilm(film, debugLines);
toc

%% Assign to optical image

oi = oiCreate;
oi = initDefaultSpectrum(oi);
oi = oiSet(oi, 'wave', wave);

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
% set(gca,'xlim',[-15 15]); set(gca,'xtick',[-15:5:15])