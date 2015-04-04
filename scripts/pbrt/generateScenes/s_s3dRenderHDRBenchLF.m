%% Render the HDR scene of the bench scene as a light field (oi)
%
% 

%%
ieInit

%% Set up the general pbrt imaging parameters 

clear curPbrt;
curPbrt = pbrtObject();
      
%make film and lens and add to pbrt object
filmDist = 45;  % 40; %36.77;
filmDiag = 12;

lensFile    = 'dgauss.50mm.dat';
focalLength = 0.050;
apertureDiameter = 16;
diffraction      = false;
chromaticAberration = false;
pinholeExitApLoc    = [2  2 -35];

% This is the microlens (pinhole) array
numPinholesW = 160;      % These 2 parameters must be even (for now)
numPinholesH = 160;
microlensMode = false;  % Feature will arrive

% Assign pinhole position to PBRT, and figure out correct cropWindow
lens = pbrtLensRealisticObject(filmDist, filmDiag, lensFile, apertureDiameter, ...
    diffraction, chromaticAberration, [], [], [], ...
    numPinholesW, numPinholesH, microlensMode);
curPbrt.camera.setLens(lens);

% Set some spatialfield of view and resolution parameters
curPbrt.camera.setResolution(300, 300);    % Low Quality mode
curPbrt.camera.setCropWindow(0, 1, 0, 1);  % Show everything

% Sampler
samples = curPbrt.sampler.removeProperty();
samples.value = 64;
curPbrt.sampler.addProperty(samples);


%% Build the bench scene specific parameters

camPos = [
  324.9416 -307.5279   51.3868
  324.2154 -306.8445   51.3121
   -0.0470    0.0591    0.9971];
curPbrt.camera.setPosition(camPos);
%curPbrt.camera.setLens(pbrtLensPinholeObject(70, 70)); 

% Backdrop Depth - for this scene, the near and far benchess
backDropDepth = -100;
foregroundDepth = -65;

% Light source
curPbrt.removeLight();
curPbrt.addLightSource(fullfile(s3dRootPath,'data','pbrtScenes','benchScene', 'sunsetEnvironmentLight.pbrt'));

% Add material file
curPbrt.addMaterial(fullfile(s3dRootPath, 'data', 'pbrtScenes', 'benchScene', 'default-mat.pbrt'));

% Remove default geometry
curPbrt.removeGeometry();

% Add geometry
curPbrt.addGeometry(fullfile(s3dRootPath, 'data', 'pbrtScenes', 'benchScene','default-geom-big-bigfloor.pbrt'));

%% If you would like to check, render the scene and depth map

% lens = pbrtLensRealisticObject(filmDist, filmDiag, specFile, ...
%     apertureDiameter, diffraction, chromaticAberration, [], []);
% curPbrt.camera.setLens(lens);
% oi = s3dRenderOIAndDepthMap(curPbrt, focalLength, 'benchLF', true);
% vcAddObject(oi); oiWindow;

% Remember to put back the other lens

%% Render the LF using the new docker container

% Reset the sampler because of some cloning limitations
sampler = pbrtSamplerObject();
tempProp = sampler.removeProperty();   %remove integer pixelsamples, modify, then add back in
tempProp.value = 64;
sampler.addProperty(tempProp);
curPbrt.setSampler(sampler);

curPbrt.camera.setResolution(1440, 1440);
dockerFlag = true;

oi = s3dRenderOIAndDepthMap(curPbrt, focalLength, 'benchLF',dockerFlag);
vcAddObject(oi); oiWindow;

%% To save do this
save('benchLF','oi','numPinholesW','numPinholesH','focalLength','curPbrt');

%% End


