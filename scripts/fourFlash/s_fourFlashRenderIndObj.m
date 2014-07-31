%% Runs PBRT and imports it in ISET for the bench scene. 
% TODO: incorporate this into 1 main script, or generalize/organize it

s_initISET

%% Light parameters

flashSeparation = 20; %mm

%% render scene with PBRT using pbrtObjects (upper right flash)
tic

% initialization
clear curPbrt;
curPbrt = pbrtObject();

% specify scene properties
matFile = fullfile(dataPath, 'twoFlashDepth', 'indObject', 'pbrt', 'graycard-mat.pbrt');
geoFile = fullfile(dataPath, 'twoFlashDepth', 'indObject', 'pbrt', 'lambertian-geom.pbrt');

% camera properties
from = [ -56.914787 -105.385544 35.0148];
to = [-56.487434 -104.481461 34.8 ];

% the vector that indicates the camera direction
camDirVector = to - from;  

% calculate a perpendicular vector (in the perpendicular plane)
perpHorVector = [camDirVector(2)  -camDirVector(1) 0];
perpHorVector = perpHorVector./sqrt(sum(perpHorVector .* perpHorVector));

%caculate a vertical vector using a cross product
perpVerVector = cross(perpHorVector, camDirVector);
perpVerVector = normvec(perpVerVector);

upVector = perpVerVector;


position = [from; to; upVector];
curPbrt.camera.setPosition(position);
lens = pbrtLensPinholeObject();
filmDistance = 140;
filmDiag = 50.9117;
curPbrt.camera.setLens(pbrtLensPinholeObject(filmDistance, filmDiag));  %TODO: may want to switch to a real lens later

% light properties
spectrum = pbrtSpectrumObject('rgb I', [1000 1000 1000]);
lightFrom = [from + perpHorVector * flashSeparation/2 + perpVerVector * flashSeparation/2];
lightTo = lightFrom + camDirVector;
coneAngle = 180;
coneDeltaAngle = 180;
lightSource = pbrtLightSpotObject('light', spectrum, coneAngle, coneDeltaAngle, lightFrom, lightTo);  %lightSpotObject(inName, inSpectrum, inConeAngle, inDeltaAngle, inFrom, inTo)

% add old parts, put in new ones
curPbrt.removeMaterial();
curPbrt.addMaterial(matFile);
curPbrt.removeGeometry();
curPbrt.addGeometry(geoFile);
curPbrt.removeLight();
curPbrt.addLightSource(lightSource);

%write file and render
tmpFileName = 'indObjUpperRightFlash';
curPbrt.writeFile(tmpFileName);
frontOi = s3dRenderOI(curPbrt, .050, tmpFileName);

toc
%% render scene with PBRT using pbrtObjects (upper Left)
tic

%light properties
lightFrom = [from + -perpHorVector * flashSeparation/2 + perpVerVector * flashSeparation/2];
lightTo = lightFrom + camDirVector;
lightSource = pbrtLightSpotObject('light', spectrum, coneAngle, coneDeltaAngle, lightFrom, lightTo);  %lightSpotObject(inName, inSpectrum, inConeAngle, inDeltaAngle, inFrom, inTo)

% add old parts, put in new ones
curPbrt.removeLight();
curPbrt.addLightSource(lightSource);

%write file and render
tmpFileName = 'indObjUpperLeftFlash';
curPbrt.writeFile(tmpFileName);
backOi = s3dRenderOI(curPbrt, .050, tmpFileName);

toc

%% render scene with PBRT using pbrtObjects (lower right)
tic

%light properties
lightFrom = [from + perpHorVector * flashSeparation/2 - perpVerVector * flashSeparation/2];
lightTo = lightFrom + camDirVector;
lightSource = pbrtLightSpotObject('light', spectrum, coneAngle, coneDeltaAngle, lightFrom, lightTo);  %lightSpotObject(inName, inSpectrum, inConeAngle, inDeltaAngle, inFrom, inTo)

% add old parts, put in new ones
curPbrt.removeLight();
curPbrt.addLightSource(lightSource);

%write file and render
tmpFileName = 'indObjLowerRightFlash';
curPbrt.writeFile(tmpFileName);
backOi = s3dRenderOI(curPbrt, .050, tmpFileName);

toc

%% render scene with PBRT using pbrtObjects (lower Left)
tic

%light properties
lightFrom = [from + -perpHorVector * flashSeparation/2 - perpVerVector * flashSeparation/2];
lightTo = lightFrom + camDirVector;
lightSource = pbrtLightSpotObject('light', spectrum, coneAngle, coneDeltaAngle, lightFrom, lightTo);  %lightSpotObject(inName, inSpectrum, inConeAngle, inDeltaAngle, inFrom, inTo)

% add old parts, put in new ones
curPbrt.removeLight();
curPbrt.addLightSource(lightSource);

%write file and render
tmpFileName = 'indObjLowerLeftFlash'
curPbrt.writeFile(tmpFileName);
backOi = s3dRenderOI(curPbrt, .050, tmpFileName);

toc

%% render depthMap with PBRT using pbrtObjects
tic

%TODO: make this into a function

%change the sampler to stratified for non-noisy depth map
samplerProp = pbrtPropertyObject();
curPbrt.sampler.setType('stratified');
curPbrt.sampler.removeProperty();
curPbrt.sampler.addProperty(pbrtPropertyObject('integer xsamples', '1'));
curPbrt.sampler.addProperty(pbrtPropertyObject('integer ysamples', '1'));
curPbrt.sampler.addProperty(pbrtPropertyObject('bool jitter', '"false"'));

%write file and render
tmpFileName = 'indObjDepthMap';
curPbrt.writeFile(tmpFileName);
groundTruthDepthMap = s3dRenderDepthMap(tmpFileName, 1);
figure; imagesc(groundTruthDepthMap);

toc

%% front flash image processing
% %load oi from file (optional)
% % vcLoadObject('opticalimage', ['50mmFront.pbrt.mat']);
% % oi = vcGetObject('oi');
% 
% % sensor processing
% sensor = s3dProcessSensor(frontOi, 0, [400 400],0, 'analog');    %low noise, auto exposure
% % sensor = s3dProcessSensor(oi, .0096, [], .03);     %high noise
% vcAddAndSelectObject('sensor',sensor); sensorImageWindow;
% 
% % image processing
% vciFlash = s3dProcessImage(sensor);
% vcAddAndSelectObject(vciFlash); vcimageWindow;
% 
% %% back flash image processing
% %load oi from file
% % vcLoadObject('opticalimage', ['50mmBack.pbrt.mat']);
% % oi = vcGetObject('oi');
% 
% % sensor processing
% frontFlashExpDur = sensorGet(sensor, 'expTime');
% sensor = s3dProcessSensor(backOi, 0, [400 400],frontFlashExpDur, 'analog');    %low noise
% % sensor = s3dProcessSensor(oi, .0096, [], .03);     %high noise
% vcAddAndSelectObject('sensor',sensor); sensorImageWindow;
% 
% % image processing
% vciFlashBack = s3dProcessImage(sensor);
% vcAddAndSelectObject(vciFlashBack); vcimageWindow;

