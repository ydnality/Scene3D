%%  Creates a pointArray using point light sources.  
% The point light sources are created using tiny spherical area lights 
% diffraction is turned off in this example!

%% initialization and directory organization
filePath = [dataPath '/validate/pbrtObject/'];
chdir(filePath);
mkdir('batchPbrtFiles');
unix('rm batchPbrtFiles/*');
unix('cp * ./batchPbrtFiles/');
chdir('batchPbrtFiles');


%% make a new pbrt object
clear curPbrt;
curPbrt = pbrtObject();
curPbrt.camera.setLens('pinhole-test.pbrt');
curPbrt.sampler.setPixelSamples(256);  % curPbrt.sampler.setPixelSamples(131070);
curPbrt.camera.setResolution(201, 201);

%% run pbrt
tmpFileName = ['pinhole'  '.pbrt'];
curPbrt.writeFile(tmpFileName);
fullfname = fullfile(filePath,'batchPbrtFiles',tmpFileName);
sceneName = tmpFileName;
scene = s3dRenderScene(fullfname, sceneName);

%% Let's see what we have
vcAddObject(scene);
sceneWindow;

% chdir('..');
%the result should be a grid of points

%% End

