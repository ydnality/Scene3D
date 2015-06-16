%% s_s3dRenderDesk.m
%
% This should render the desk scene for flash/noflash cases.
%
% We use the basic format in bench and metronome files to do the rendering,
% but for the desk/default.pbrt.

%%
ieInit

%% PBRT will run the PBRT script
pbrtDir = fullfile(s3dRootPath, 'data', 'pbrtScenes','desk');
if ~exist(pbrtDir,'dir'), error('Desk pbrt directory is missing'); end

% Governs the number of samples per pixel. A minimum is 32.
% For the scarlet one we ran it with 4096.
nSamples = 1024; %512;

%%

%initialization
clear curPbrt;
curPbrt = pbrtObject();

% Specify properties files (materials and geometry)
matFile = fullfile(pbrtDir, 'default-mat.pbrt');
geoFile = fullfile(pbrtDir, 'default-geom.pbrt');


%% Standard camera properties

% This is from the default.pbrt file for desk.
LookAt  = ...
[ -127.5649  -89.6985   80.9301;
   -81.3417  -63.8818   72.3964;
     0 0 1];
curPbrt.camera.setPosition(LookAt);

% light properties
spectrum = pbrtSpectrumObject('rgb I', [1000 1000 1000]);
% lightFrom = [  -56.914787 -105.385544 35.0148];  % Position of source
% lightTo =   [-56.487434 -104.481461 34.8  ];       % Direction of principal ray
coneAngle      = 180;    % Angle of rays from light source
coneDeltaAngle = 180;    % Drop off of illumination???

lightSource = pbrtLightSpotObject('light', spectrum, coneAngle, coneDeltaAngle, LookAt(1,:), LookAt(2,:));  


%% specify the wide field camera 
    
% From the default.pbrt file
curPbrt.camera.setLens(pbrtLensRealisticObject());   
curPbrt.camera.lens.filmDistance = 52; % 70;  % 133.33;
curPbrt.camera.lens.filmDiag = 43.3;
curPbrt.camera.lens.specFile = 'dgauss.50mm.dat';
curPbrt.camera.lens.apertureDiameter = 3; % in mm
curPbrt.camera.lens.curveRadius = 0;       % Experimental
curPbrt.camera.setResolution(450, 300);
% curPbrt.camera.setResolution(64, 32);

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
curPbrt.sampler.addProperty(pbrtPropertyObject('integer pixelsamples', nSamples));

%% Render the oi
oiName = 'desk';
dockerFlag = true;
oi = s3dRenderOIAndDepthMap(curPbrt, oiName, dockerFlag);

vcAddAndSelectObject(oi);
oiWindow;

%% END

