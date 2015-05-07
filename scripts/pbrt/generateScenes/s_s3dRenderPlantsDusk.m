%% Render the Plants Dusk scene 
%  
%  Read the comments to find the key place where those parameters are set.
%
%  There is also a lightfield camera version of this rendering.
%
% AL VISTASOFT, 2014

%%
ieInit

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

% from = [28000 1800 1500 ];  %large scale version
% to =   [69000 65000 1300];

vector = to - from;
normVector = vector./sqrt(sum(vector .* vector));

from = from - normVector * 2;
to = to - normVector * 2;

% from = [33 -50 40];   %standing way back
% to = [50 50 1];

position = [from; to; 0 0 1];
curPbrt.camera.setPosition(position);
curPbrt.camera.lens.filmDistance = 30;  %increase the FOV of the pinhole camera
curPbrt.camera.lens.filmDistance = 50;  %increase the FOV of the pinhole camera


% curPbrt.camera.setLens(pbrtLensRealisticObject());   
% curPbrt.camera.lens.filmDistance = 90; % 70;  % 133.33;
% curPbrt.camera.lens.filmDiag = 70;
% curPbrt.camera.lens.specFile = '2ElLensSmallAp.dat';
% curPbrt.camera.lens.apertureDiameter = .01; % in mm
% curPbrt.camera.lens.curveRadius = 0;       % Experimental
% curPbrt.camera.setResolution(450, 300);


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
%oi = s3dRenderOI(curPbrt, oiName, dockerFlag);


toc

vcAddAndSelectObject(oi);
oiWindow;
