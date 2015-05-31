%% Runs PBRT and imports it in ISET for the bench scene. 
% TODO: incorporate this into 1 main script, or generalize/organize it

ieInit
%% render scene with PBRT using pbrtObjects (front flash)
tic

%initialization
clear curPbrt;
curPbrt = pbrtObject();

%specify scene properties
%matFile = fullfile(dataPath, 'twoFlashDepth', 'indObject', 'pbrt', 'graycard-mat.pbrt');
matFile = fullfile(dataPath, 'twoFlashDepth', 'indObject', 'pbrt', 'default-mat.pbrt');
%matFile = fullfile(dataPath, 'twoFlashDepth', 'indObject', 'pbrt', 'mixed-mat.pbrt');
%matFile = fullfile(dataPath, 'twoFlashDepth', 'indObject', 'pbrt', 'whiteWallsFloor-mat.pbrt');

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
filmDiag = 50.9117;   %36 width and height
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
frontOi = s3dRenderOIAndDepthMap(curPbrt); %, .050);
%frontOi = oiSet(frontOi,'optics focal length',0.050);  % I set these variables by hand just to avoid the zero condition
frontOi = oiSet(frontOi,'optics focal length',0.140);  % I set these variables by hand just to avoid the zero condition
frontOi = oiSet(frontOi,'optics f number', 2);    

toc

%% front flash image processing
%load oi from file (optional)
% vcLoadObject('opticalimage', ['50mmFront.pbrt.mat']);
% oi = vcGetObject('oi');

% % sensor processing
noiseFlag = 0;   %0 is off
sensor = s3dProcessSensor(frontOi, 0, [400 400],0, 'analog', noiseFlag);    %low noise, auto exposure
% sensor = s3dProcessSensor(oi, .0096, [], .03);     %high noise
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

% image processing
vciFlash = s3dProcessImage(sensor);
vcAddObject(vciFlash); ipWindow;


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
%position = [from; to; 0 0 1];
%curPbrt.camera.setPosition(position);

coneAngle = 180;
coneDeltaAngle = 180;
lightSource = pbrtLightSpotObject('light', spectrum, coneAngle, coneDeltaAngle, lightBackFrom, lightBackTo);  %lightSpotObject(inName, inSpectrum, inConeAngle, inDeltaAngle, inFrom, inTo)

% add old parts, put in new ones
curPbrt.removeLight();
curPbrt.addLightSource(lightSource);

%write file and render
backOi = s3dRenderOI(curPbrt);
backOi = oiSet(backOi,'optics focal length',0.140);  % I set these variables by hand just to avoid the zero condition
backOi = oiSet(backOi,'optics f number', 2);    

toc


%% back flash image processing
%load oi from file
% vcLoadObject('opticalimage', ['50mmBack.pbrt.mat']);
% oi = vcGetObject('oi');

% sensor processing
frontFlashExpDur = sensorGet(sensor, 'expTime');
noiseFlag = 0;
sensor = s3dProcessSensor(backOi, 0, [400 400],frontFlashExpDur, 'analog', noiseFlag);    %low noise
% sensor = s3dProcessSensor(oi, .0096, [], .03);     %high noise
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

% image processing
vciFlashBack = s3dProcessImage(sensor);
vcAddAndSelectObject(vciFlashBack); ipWindow;

% 
% frontFlashExpDur = sensorGet(sensor, 'expTime');
% 
% ss = oiGet(backOi,'sample spacing','m');
% sensor = sensorCreate;
% sensor = sensorSet(sensor,'pixel size same fill factor',ss(1));
% sensor = sensorSet(sensor,'size',oiGet(backOi,'size'));
% sensor = sensorSet(sensor,'exp time',frontFlashExpDur);
% %sensor = sensorSet(sensor,'exp time',0.25386);
% %sensor = sensorSet(sensor, 'autoexposure', 1);
% 
% % pixel = sensorGet(sensor, 'pixel')
% %sensor = sensorSet(sensor, 'prnu level', .001);
% sensor = sensorSet(sensor, 'dsnu level', .000);
% 
% % Describe
% sensorGet(sensor,'pixel size','um')
% sensorGet(sensor,'size')
% sensorGet(sensor,'fov',[],backOi)
% 
% % Compute the sensor response
% sensor = sensorCompute(sensor,backOi);
% vcAddObject(sensor); sensorWindow('scale',1);
% 
% % Interpolate the color filter data to produce a full sensor
% vciFlashBack = ipCreate;
% vciFlashBack = ipCompute(vciFlashBack,sensor);
% vcAddObject(vciFlashBack); ipWindow;