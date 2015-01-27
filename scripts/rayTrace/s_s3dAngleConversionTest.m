%% Testing of the spherical angle conversion function in the rayC class
%
% This script renders a PSF and in the process tests the spherical angle
% conversion function in the rayC class
%
% AL Vistalab 2014

%% Initialize Iset
s_initISET

%% Declare ray-trace type

rtType = 'realistic';  %ideal/realistic
debugLines = 50;
% (lowerLeftx,LowerLefty,uppRightx,upperRighty)
% Percentage of rectangular lens, 0 is middle
subSection = [];   % Whole thing
% subSection = [-.25 -.25 .25 .25];   % Not working
%subSection = [0 0 .25 .25];   

%% Declare point sources
% declare point sources in world space.  The camera is usually at [0 0 0],
% and pointing towards -z.  We are using a right-handed coordinate system.

% [XGrid YGrid] = meshgrid(-4000:1000:4000,-4000:1000:4000);
% [XGrid YGrid] = meshgrid(-2000:1000:2000,-2000:1000:2000);

% Use the psCreate thing ...
%pointSources = [0 10 -2000];
pointSources = [0 0 -100];

%% Declare camera properties

% Build a sensor (film) object
% Position, size,  wave, waveConversion, resolution
% film = filmC([0 0 51.2821	],[.2/sqrt(2) .2/sqrt(2)], 400:10:700, [(400:10:700)' (1:31)'], [50 50 31]);

% Declare film
% filmPosition = [0 0 51.2821	];  % Good for 2Elens
% filmPosition = [0 0 37.4];        % Good for dgauss.50mm.  True focal about 37.3mm
filmPosition = [0 0 107]; 


% filmSize = [.2/sqrt(2) .2/sqrt(2)];
%filmDiag = 1;  % Millimeters
filmDiag = 3;  % Millimeters
filmSize = [filmDiag/sqrt(2) filmDiag/sqrt(2)];

wave = 400:50:600;
resolution =  [300 300 length(wave)];
film = filmC('position', filmPosition, 'size', filmSize, 'wave', wave, 'resolution', resolution);

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

% Make a function that goes gets a lens from a file, say
% lens = lensC('filename',fname);
%  That should go and see if there is a file called fname.
%
lensFile = fullfile(s3dRootPath, 'data', 'lens', '2ElLens.mat');
import = load(lensFile,'lens');
thickLens = import.lens;
thickLens.apertureMiddleD = 10;

lensFile = fullfile(s3dRootPath, 'data', 'lens', 'dgauss.50mm.mat');
import = load(lensFile,'lens');
multiLens = import.lens;

% lens = multiLens;
lens = multiLens;
lens.set('wave', wave);
% Matrix of n for each surface element.
% Apertures are 0
n = lens.get('nArray');
% Plot it?  Make an image?  Have a swell time

%% calculate the origin and direction of the rays
curInd = 1;
disp('-----trace source to lens-----');
tic
rays = lens.rtSourceToEntrance(pointSources(curInd, :), false, jitterFlag, rtType, subSection)
toc

%  look at angles and analyze them
%sphereAngles = rays.get('sphericalAngles');

% figure; hist(sphereAngles(:,1));  %should be uniform
%should have more at edges (more data points at perimeter) ... approximately linear
% figure; hist(sphereAngles(:,2));  

%% use the projection form of angles and plot phase space
projAngles = rays.get('projectedAngles');
% hist(projAngles(:,1)); 
% hist(projAngles(:,2)); 

% plot phase space
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

%% 