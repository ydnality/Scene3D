%% Runs PBRT and imports it in ISET for the bench scene. 
% TODO: incorporate this into 1 main script, or generalize/organize it

s_initISET
%% render scene with PBRT using pbrtObjects (front flash)
tic

%initialization
clear curPbrt;
curPbrt = pbrtObject();

%specify scene properties
%matFile = fullfile(dataPath, 'twoFlashDepth', 'indObject', 'pbrt', 'graycard-mat.pbrt');
%matFile = fullfile(dataPath, 'twoFlashDepth', 'indObject', 'pbrt', 'default-mat.pbrt');
matFile = fullfile(dataPath, 'twoFlashDepth', 'indObject', 'pbrt', 'mixed-mat.pbrt');

%geoFile = fullfile(dataPath, 'twoFlashDepth', 'indObject', 'pbrt', 'lambertian-geom.pbrt');
geoFile = fullfile(dataPath, 'twoFlashDepth', 'indObject', 'pbrt', 'default-geom.pbrt');


%light properties
spectrum = pbrtSpectrumObject('rgb I', [1000 1000 1000]);
lightFrom = [  -56.914787 -105.385544 35.0148];
lightTo = [-56.487434 -104.481461 34.8  ];
coneAngle = 180;
coneDeltaAngle = 180;
lightSource = pbrtLightSpotObject('light', spectrum, coneAngle, coneDeltaAngle, lightFrom, lightTo);  %lightSpotObject(inName, inSpectrum, inConeAngle, inDeltaAngle, inFrom, inTo)

%camera properties
from = [ -56.914787 -105.385544 35.0148];
to = [-56.487434 -104.481461 34.8 ];
position = [from; to; 0 0 1];
curPbrt.camera.setPosition(position);
lens = pbrtLensPinholeObject();
filmDistance = 140;
filmDiag = 50.9117;
curPbrt.camera.setLens(pbrtLensPinholeObject(filmDistance, filmDiag));  %TODO: may want to switch to a real lens later
curPbrt.camera.setResolution(300, 300);

% add old parts, put in new ones
curPbrt.removeMaterial();
curPbrt.addMaterial(matFile);
curPbrt.removeGeometry();
curPbrt.addGeometry(geoFile);
curPbrt.removeLight();
curPbrt.addLightSource(lightSource);

% set sampler
curPbrt.sampler.removeProperty();
curPbrt.sampler.addProperty(pbrtPropertyObject('integer pixelsamples', 256));

%write file and render
frontOi = s3dRenderOI(curPbrt, .050);

toc
%% render scene with PBRT using pbrtObjects (back flash)
tic

%** must run first part of last section for now.  Cloning while rendering
%is broken.

%light properties
spectrum = pbrtSpectrumObject('rgb I', [1000 1000 1000]);
lightBackFrom = [ -77.8060 -149.5817   45.5153];
lightBackTo = [-77.3786 -148.6776   45.3005 ];
from = lightBackFrom;
to = lightBackTo;
position = [from; to; 0 0 1];
curPbrt.camera.setPosition(position);

coneAngle = 180;
coneDeltaAngle = 180;
lightSource = pbrtLightSpotObject('light', spectrum, coneAngle, coneDeltaAngle, lightBackFrom, lightBackTo);  %lightSpotObject(inName, inSpectrum, inConeAngle, inDeltaAngle, inFrom, inTo)

% add old parts, put in new ones
curPbrt.removeLight();
curPbrt.addLightSource(lightSource);

%write file and render
backOi = s3dRenderOI(curPbrt, .050);

toc
%% render depthMap with PBRT using pbrtObjects
tic

%change the sampler to stratified for non-noisy depth map
samplerProp = pbrtPropertyObject();
curPbrt.sampler.setType('stratified');
curPbrt.sampler.removeProperty();
curPbrt.sampler.addProperty(pbrtPropertyObject('integer xsamples', '1'));
curPbrt.sampler.addProperty(pbrtPropertyObject('integer ysamples', '1'));
curPbrt.sampler.addProperty(pbrtPropertyObject('bool jitter', '"false"'));

%write file and render
tmpFileName = ['deleteMe'  '.pbrt'];
curPbrt.writeFile(tmpFileName);
groundTruthDepthMap = s3dRenderDepthMap(tmpFileName, 1);
figure; imagesc(groundTruthDepthMap);

toc
%% front flash image processing
%load oi from file (optional)
% vcLoadObject('opticalimage', ['50mmFront.pbrt.mat']);
% oi = vcGetObject('oi');

% sensor processing
sensor = s3dProcessSensor(frontOi, 0, [400 400],0, 'analog');    %low noise, auto exposure
% sensor = s3dProcessSensor(oi, .0096, [], .03);     %high noise
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

% image processing
vciFlash = s3dProcessImage(sensor);
vcAddAndSelectObject(vciFlash); vcimageWindow;

%% back flash image processing
%load oi from file
% vcLoadObject('opticalimage', ['50mmBack.pbrt.mat']);
% oi = vcGetObject('oi');

% sensor processing
frontFlashExpDur = sensorGet(sensor, 'expTime');
sensor = s3dProcessSensor(backOi, 0, [400 400],frontFlashExpDur, 'analog');    %low noise
% sensor = s3dProcessSensor(oi, .0096, [], .03);     %high noise
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

% image processing
vciFlashBack = s3dProcessImage(sensor);
vcAddAndSelectObject(vciFlashBack); vcimageWindow;

%%  generate and read and output depth map (obsolete)

% % ** make sure the rendering file has a small initial aperture, and only 1
% % sample per pixel!
% 
% %currently broken - works - but we need to change the sampler to
% %"stratified" with 1 sample in each direction
% 
% % chdir(fullfile(datapath, 'twoFlashDepth', 'depthTargetDepths'));
% chdir(fullfile(datapath, 'twoFlashDepth', 'indObject', 'pbrt'));
% 
% depthMap = s3dRenderDepthMap('50mmDepthMap.pbrt', 1);
% figure; imagesc(depthMap);

%%  generate and read and output depth map (obsolete)
% ** make sure the rendering file has a small initial aperture, and only 1
% sample per pixel!

%works - but we need to change the sampler to
%"stratified" with 1 sample in each direction

% % chdir(fullfile(datapath, 'twoFlashDepth', 'depthTargetDepths'));
% chdir(fullfile(datapath, 'twoFlashDepth', 'indObject', 'pbrt'));
% 
% depthMap = s3dRenderDepthMap('50mmDepthMap.pbrt', 1);
% figure; imagesc(depthMap);

%% render scene with PBRT using classic pbrt files(obsolete)

% %new batch render code
% tic
% 
% %define template file
% % templateFile = '2FlashDepthTemplate';  %(no extension on here)
% templateFile = 'templateDepthTargetDepths';  %(no extension on here)
% 
% 
% clear('templateFileArray');
% templateFileArray{1} = batchFileClass(templateFile, '.pbrt');
% % templateFileArray{2} = batchFileClass(templateFile, '-idealLens.pbrt');
% % templateFileArray{3} = batchFileClass(templateFile, '-realisticLens.pbrt');
% % templateFileArray{4} = batchFileClass(templateFile, '-frontFlash.pbrt');
% % templateFileArray{5} = batchFileClass(templateFile, '-backFlash.pbrt');
% %define conditions file
% 
% % conditionsFile = 'conditionsScratch.txt';
% conditionsFile = 'conditions10mSep.txt';
% 
% %change into correct directory
% %chdir(twoFlashDepthPath);
% chdir([dataPath '/twoFlashDepth/depthTargetDepths/']);
% renderFileArray = s3dRenderBatch(templateFileArray, conditionsFile);
% 
% toc
