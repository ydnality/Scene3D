%% Runs PBRT and imports it in ISET for the bench scene. 

%% render scene with PBRT
sceneName = 'indObjAmbientOnly';
% oi = s3dRenderOI('indestructibleObject/simpleReflectance-downAmbient.pbrt', .050);
% oi = s3dRenderOI('indestructibleObject/graycard-down.pbrt', .050);
% oi = s3dRenderOI('desk/paintWhiteFlash.pbrt', .050);
oi = s3dRenderOI('desk/graycardWhiteFlash.pbrt', .050);
% oi = s3dRenderOI('indestructibleObject/paintReflectance-downWhiteFlash.pbrt', .050);
% oi = s3dRenderOI('indestructibleObject/simpleReflectance-downAmbientOnly.pbrt', .050);

oi = oiSet(oi, 'name', sceneName);
% strip the file name from the path and assign that as the name of the
% object  ... vcSaveObject(oi,);
vcAddAndSelectObject(oi);
oiWindow;

%% camera processing (sensor and image) - no flash image
% run this section for the no flash image

% sensor processing
sensor = s3dProcessSensor(oi, 0, [], 0, 'analog');    %low noise, auto exposure
% sensor = s3dProcessSensor(oi, .0096, [], .03);     %high noise
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

% image processing
image = s3dProcessImage(sensor);
vcAddAndSelectObject(image); vcimageWindow;

%% camera processing (sensor and image) -  flash image
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
depthMap = s3dRenderDepthMap('indestructibleObject/simpleReflectance-downAmbient.pbrt', 11);
figure; imagesc(depthMap);