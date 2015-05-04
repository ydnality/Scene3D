%% Render the Metronome scene with a standard camera or wide field
%
%  With some parameter settings below, we also can render the scene with a
%  very wide field of view and with a curved sensor. 
%  
%  Read the comments to find the key place where those parameters are set.
%
%  There is also a lightfield camera version of this rendering.
%
% AL VISTASOFT, 2014

%%
ieInit

%% Rendering parameters
%wideField    = false;


%%
tic

%initialization
clear curPbrt;
curPbrt = pbrtObject();

% Specify properties files (materials and geometry)
matFile = fullfile(dataPath, 'pbrtScenes', 'plantsDusk', 'ecosys-texture.pbrt');

geoFile1 = fullfile(dataPath, 'pbrtScenes', 'plantsDusk', 'ecosys-apple.pbrt');
geoFile2 = fullfile(dataPath, 'pbrtScenes', 'plantsDusk', 'ecosys-bush.pbrt');
geoFile3 = fullfile(dataPath, 'pbrtScenes', 'plantsDusk', 'ecosys-dandelion.pbrt');
geoFile4 = fullfile(dataPath, 'pbrtScenes', 'plantsDusk', 'ecosys-grass.pbrt');
geoFile5 = fullfile(dataPath, 'pbrtScenes', 'plantsDusk', 'ecosys-urtica.pbrt');
geoFile6 = fullfile(dataPath, 'pbrtScenes', 'plantsDusk', 'ecosys-yellowflower.pbrt');
geoFile7 = fullfile(dataPath, 'pbrtScenes', 'plantsDusk', 'ecosys-terrain.pbrt');


%% Standard camera properties
from = [28 1.8 1.5 ];
to =   [69 65 1.3];
position = [from; to; 0 0 1];
curPbrt.camera.setPosition(position);
curPbrt.camera.lens.filmDistance = 30;  %increase the FOV of the pinhole camera

%% add old parts, put in new ones

% Candidate for a function
curPbrt.removeMaterial();
curPbrt.addMaterial(matFile);
curPbrt.removeGeometry();
curPbrt.addGeometry(geoFile1);
curPbrt.addGeometry(geoFile2);
curPbrt.addGeometry(geoFile3);
curPbrt.addGeometry(geoFile4);
curPbrt.addGeometry(geoFile5);
curPbrt.addGeometry(geoFile6);
curPbrt.addGeometry(geoFile7);
curPbrt.removeLight();
%curPbrt.addLightSource(lightSource);

includeFile = fullfile(dataPath, 'pbrtScenes', 'plantsDusk',  'plantsDusk-include.pbrt');
curPbrt.addInclude(includeFile);

% set sampler
curPbrt.sampler.removeProperty();
nSamples = 512; %512;
curPbrt.sampler.addProperty(pbrtPropertyObject('integer pixelsamples', nSamples));

%set volume integrator to default
curPbrt.setVolumeIntegrator(pbrtVolumeIntegratorObject());

%% Render the oi
dockerFlag = true;
oiName = 'plantsDusk';
oi = s3dRenderOIAndDepthMap(curPbrt, oiName, dockerFlag);

toc

vcAddAndSelectObject(oi);
oiWindow;
