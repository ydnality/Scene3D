%% Runs PBRT and imports it in ISET for the bench scene. 

% TODO: incorporate this into 1 main script, or generalize/organize it

%% render scene with PBRT

%old stuff
% sceneName = 'indObjAmbientOnly';
% oi = s3dRenderScene('indestructibleObject/graycard-down.pbrt', .050);
% 
% oi = oiSet(oi, 'name', sceneName);
% % strip the file name from the path and assign that as the name of the
% % object  ... vcSaveObject(oi,);
% vcAddAndSelectObject(oi);
% oiWindow;


%new batch render code
tic

%define template file
% templateFile = '2FlashDepthTemplate';  %(no extension on here)
templateFile = 'templateDepthTargetDepths';  %(no extension on here)


clear('templateFileArray');
templateFileArray{1} = batchFileClass(templateFile, '.pbrt');
% templateFileArray{2} = batchFileClass(templateFile, '-idealLens.pbrt');
% templateFileArray{3} = batchFileClass(templateFile, '-realisticLens.pbrt');
% templateFileArray{4} = batchFileClass(templateFile, '-frontFlash.pbrt');
% templateFileArray{5} = batchFileClass(templateFile, '-backFlash.pbrt');
%define conditions file

% conditionsFile = 'conditionsScratch.txt';
conditionsFile = 'conditions10mSep.txt';

%change into correct directory
%chdir(twoFlashDepthPath);
chdir([dataPath '/twoFlashDepth/depthTargetDepths/']);
renderFileArray = s3dRenderBatch(templateFileArray, conditionsFile);

toc

%% front flash image processing

%load oi from file
% chdir(TwoFlashDepthPath);
vcLoadObject('opticalimage', ['50mmFront.pbrt.mat']);
oi = vcGetObject('oi');


% sensor processing
% sensor = s3dProcessSensor(oi, 0, [],.356, 'analog');    %low noise, auto exposure
sensor = s3dProcessSensor(oi, 0, [],.105, 'analog');    %low noise, auto exposure

% sensor = s3dProcessSensor(oi, .0096, [], .03);     %high noise
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

% image processing
image = s3dProcessImage(sensor);
vcAddAndSelectObject(image); vcimageWindow;


%% back flash image processing
%load oi from file
% chdir(TwoFlashDepthPath);
vcLoadObject('opticalimage', ['50mmBack.pbrt.mat']);
oi = vcGetObject('oi');


% sensor processing
% sensor = s3dProcessSensor(oi, 0, [],.356, 'analog');    %low noise, auto exposure
sensor = s3dProcessSensor(oi, 0, [],.105, 'analog');    %low noise, auto exposure
% sensor = s3dProcessSensor(oi, .0096, [], .03);     %high noise
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

% image processing
image = s3dProcessImage(sensor);
vcAddAndSelectObject(image); vcimageWindow;


%% camera processing (sensor and image) - no flash image (legacy stuff)
% run this section for the no flash image



% sensor processing
sensor = s3dProcessSensor(oi, 0, [], 0, 'analog');    %low noise, auto exposure
% sensor = s3dProcessSensor(oi, .0096, [], .03);     %high noise
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

% image processing
image = s3dProcessImage(sensor);
vcAddAndSelectObject(image); vcimageWindow;



% camera processing (sensor and image) -  flash image
% run this section for the flash image

% Uncomment this section if you wish to load in a saved oi
% load('compPhotography/reflectanceRecovery/indObjSimpleRadianceFrontFlashOi.mat'); 
% oi = opticalimage;


% sensor processing
sensor = s3dProcessSensor(oi, 0, [],  .08, 'analog');  %we might need to play with exposure
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

%image processing
[image, transformMatrices] = s3dProcessImage(sensor, []);
vcAddAndSelectObject(image); vcimageWindow;

%%  generate and read and output depth map
% ** make sure the rendering file has a small initial aperture, and only 1
% sample per pixel!
chdir(TwoFlashDepthPath);
depthMap = s3dRenderDepthMap('50mmDepthMap.pbrt', 11);
figure; imagesc(depthMap);
