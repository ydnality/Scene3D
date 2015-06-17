%% Runs PBRT and to create an ISET optical image for the metronome scene
%
% This is sometimes called indestructible object, for some reason.
%
% See also: s_s3dRenderHDRBenchLF.m
%
% AL Vistasoft 2015

%%
ieInit

%% render scene with PBRT using pbrtObjects (front flash)
tic

%pbrt object initialization
clear curPbrt;
curPbrt = pbrtObject();

% Specify properties files (materials and geometry)
matFile = fullfile(dataPath, 'twoFlashDepth', 'indObject', 'pbrt', 'default-mat.pbrt');
geoFile = fullfile(dataPath, 'twoFlashDepth', 'indObject', 'pbrt', 'default-geom.pbrt');

% light properties
spectrum = pbrtSpectrumObject('rgb I', [1000 1000 1000]);
lightFrom = [ -56.914787 -105.385544 35.0148];  % Position of source
lightTo =   [-56.487434  -104.481461 34.8  ];   % Direction of principal ray
coneAngle      = 180;    % Angle of rays from light source
coneDeltaAngle = 180;    % Drop off of illumination???
lightSource = pbrtLightSpotObject('light', spectrum, coneAngle, coneDeltaAngle, lightFrom, lightTo);  

% lightSpotObject(inName, inSpectrum, inConeAngle, inDeltaAngle, inFrom, inTo)

%% camera properties
from = [ -56.914787 -105.385544 35.0148];
to = [-56.487434 -104.481461 34.8 ];
position = [from; to; 0 0 1];
curPbrt.camera.setPosition(position);

%% add old parts, put in new ones
% Candidate for a function
curPbrt.removeMaterial();
curPbrt.addMaterial(matFile);
curPbrt.removeGeometry();
curPbrt.addGeometry(geoFile);
curPbrt.removeLight();
curPbrt.addLightSource(lightSource);

% set sampler
curPbrt.sampler.removeProperty();
nSamples = 64; %512;
curPbrt.sampler.addProperty(pbrtPropertyObject('integer pixelsamples', nSamples));

%% setup the lightfield camera properties

numPinholesW = 160;  %these 2 parameters must be even (for now)
numPinholesH = 160;
rows = 1440; cols = 1440;

filmDist = 72; %40; %36.77; mm from the back of the lens
filmDiag = 30; % Sensor diagonal size

% Lens propertiess
specFile = 'dgauss.50mm.dat';  % The lens file % specFile = '2ElLens.dat';
apertureDiameter = 16;         % Units?
diffraction = false;
chromaticAberration =false;

%assign pinhole position to PBRT, and figure out correct cropWindow
lens = pbrtLensRealisticObject(filmDist, filmDiag, ...
    specFile, apertureDiameter, ...
    diffraction, chromaticAberration, [], [], [], ...
    numPinholesW, numPinholesH);
curPbrt.camera.setLens(lens);

curPbrt.camera.setResolution(rows, cols);

%% To use a pinhole, without the microlens array set the lens this way
%  The script s_s3dRenderMetronome.m also has examples

%lens         = pbrtLensPinholeObject();
%filmDistance = 140;
%filmDiag     = 50.9117;
%TODO: may want to switch to a real lens later
%curPbrt.camera.setLens(pbrtLensPinholeObject(filmDistance, filmDiag)); 
%curPbrt.camera.setResolution(300, 300);

%% write file and render
focalLength = 0.050;   % In meters
oiName = [];
dockerFlag = true;

oi = s3dRenderOIAndDepthMap(curPbrt,oiName,dockerFlag);
vcAddObject(oi); oiWindow;

%% Save all the key variables that enable running the previous cell
save('metronomeLF','oi','numPinholesW','numPinholesH','focalLength','curPbrt');

%% END
