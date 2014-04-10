%%  Creates a pointArray using point light sources.  
% The point light sources are created using tiny spherical area lights 
% diffraction is turned off in this example!

%% make a new pbrt object
clear curPbrt;
curPbrt = pbrtObject();
% curPbrt.camera.setLens('pinhole-test.pbrt');
curPbrt.sampler.removeProperty();
curPbrt.sampler.addProperty(pbrtPropertyObject('integer pixelsamples', 256));
% curPbrt.sampler.setPixelSamples(256);  % curPbrt.sampler.setPixelSamples(131070);
curPbrt.camera.setResolution(201, 201);

%% run pbrt
scene = s3dRenderSceneAndDepthMap(curPbrt, 'pinhole');
vcAddAndSelectObject(scene);
sceneWindow;


%old stuff
% 
% curPbrt.writeFile(fullfname);
% scene = s3dRenderScene(fullfname, sceneName);
% 
% %% depth map - consider putting this into it's own function
% 
% %inputs - pbrt object, output, rendered scene, and rendered depth map
% %the lens must be a pinhole camera, or a lens with a pinhole aperture
% curPbrt.sampler.setType('stratified');
% curPbrt.sampler.removeProperty();
% curPbrt.sampler.addProperty(pbrtPropertyObject('bool jitter', '"false"'));
% curPbrt.sampler.addProperty(pbrtPropertyObject('integer xsamples', 1));
% curPbrt.sampler.addProperty(pbrtPropertyObject('integer ysamples', 1));
% 
% tmpFileName = ['pinholeDepthMap'  '.pbrt'];
% curPbrt.writeFile(tmpFileName);
% fullfnameDM = fullfile(filePath,'batchPbrtFiles',tmpFileName);
% depthMap = s3dRenderDepthMap(fullfnameDM, 1);  %this is in mm
% scene = sceneSet(scene, 'depthMap', depthMap);
% 
% %% Let's see what we have
% vcAddAndSelectObject(scene);
% sceneWindow;

% chdir('..');
%the result should be a grid of points

%% End

