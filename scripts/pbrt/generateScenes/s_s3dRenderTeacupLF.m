%% Runs PBRT and to create an ISET optical image for the metronome scene
%
% This is sometimes called indestructible object, for some reason.
%
% See also: s_s3dRenderHDRBenchLF.m, s_*MetronomeLF.m
% 
%  Not running yet -- docker container runs, but oi is black.
%
% BW Vistasoft 2015

%%
ieInit

%% render scene with PBRT using pbrtObjects (front flash)

clear curPbrt;
curPbrt = pbrtObject();

% Specify properties files (materials and geometry)
% matFile = fullfile(dataPath, 'twoFlashDepth', 'indObject', 'pbrt', 'default-mat.pbrt');
% geoFile = fullfile(dataPath, 'twoFlashDepth', 'indObject', 'pbrt', 'default-geom.pbrt');
matFile = fullfile(dataPath, 'pbrtScenes', 'teacup', 'teacup-mat.pbrt');
geoFile = fullfile(dataPath, 'pbrtScenes', 'teacup', 'teacup-geom.pbrt');

% light properties - DOn't understand all of these, but in particular where
% the direction came from.
spectrum = pbrtSpectrumObject('rgb I', [1000 1000 1000]);
lightFrom = [ 220.784576 -263.762909 12.850847];
lightTo =   [ 128.038239 -200.766434 -18.473597];
% lightFrom = [-56.914787 -105.385544 35.0148];  % Position of source
% lightTo =   [-56.487434 -104.481461 34.8  ];       % Direction of principal ray
coneAngle      = 180;    % Angle of rays from light source
coneDeltaAngle = 180;    % Drop off of illumination???
lightSource = pbrtLightSpotObject('light', spectrum, coneAngle, coneDeltaAngle, lightFrom, lightTo);  

%% camera properties

% This is the camera position, taken from the teacup.pbrt file.  I think
% that must be the 'LookAt' parameter.  Notice that they are the same as
% the light positions, above.
from = [ 220.784576 -263.762909 12.850847];
to =   [ 128.038239 -200.766434 -18.473597];
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
nSamples = 32; %512;
curPbrt.sampler.addProperty(pbrtPropertyObject('integer pixelsamples', nSamples));

%% setup the lightfield camera properties

numPinholesW = 160/2;  %these 2 parameters must be even (for now)
numPinholesH = 160/2;
rows = 1440/2; cols = 1440/2;

filmDist = 75;   % teacup.pbrt - 40; %36.77; mm from the back of the lens
filmDiag = 43.3; % teacup.pbrt - Sensor diagonal size

% Lens propertiess
specFile = 'dgauss.50mm.dat';  % teacup.pbrt - The lens file % specFile = '2ElLens.dat';
apertureDiameter = 16;         % teacup.pbrt
diffraction         = false;
chromaticAberration = false;

%assign pinhole position to PBRT, and figure out correct cropWindow
lens = pbrtLensRealisticObject(filmDist, filmDiag, ...
    specFile, apertureDiameter, ...
    diffraction, chromaticAberration, [], [], [], ...
    numPinholesW, numPinholesH);
curPbrt.camera.setLens(lens);

curPbrt.camera.setResolution(rows, cols);

%% write file and render
focalLength = 0.050;   % In meters
oiName = [];
dockerFlag = true;

oi = s3dRenderOIAndDepthMap(curPbrt,focalLength,oiName,dockerFlag);
vcAddObject(oi); oiWindow;

%% Save all the key variables that enable running the previous cell
save('teacup','oi','numPinholesW','numPinholesH','focalLength','curPbrt');

%% END
