%% Render Scene Radiance Using pbrtObjects
%This particular scene is a camera within a sphere.  
%In retrospect, this is a very boring scene because the depth will just be
%the same everywhere, since we are imaging a sphere!

s_initISET
%% 
clear curPbrt;
curPbrt = pbrtObject();

%% camera 

% camera position
newCamPos =    [0  0 0;
    0   0 -1;
    0 1.00000 0];
curPbrt.camera.setPosition(newCamPos);
curPbrt.camera.lens.filmDistance = 133.33;
curPbrt.camera.lens.filmDiag = 70;
% curPbrt.camera.setResolution(100, 100);    %LQ mode

% uncomment to use a 2 element lens instead of a pinhole
% curPbrt.camera.setLens(fullfile(s3dRootPath, 'data', 'lens', '2ElLens50mm.pbrt'));

% sampler
sampler = curPbrt.sampler.removeProperty();
sampler.value = 128;
curPbrt.sampler.addProperty(sampler);

%% Light Sources

% simple spot light object
curPbrt.removeLight();
lightFront = pbrtLightSpotObject('lightFront', [], [], [], [0 0 0], [0 0 -1]);
curPbrt.addLightSource(lightFront);

%% Scene Materials

% add a new material
matRGB= [1 1 1];
newMaterial = pbrtMaterialObject('grayMat', 'matte', pbrtPropertyObject('color Kd', matRGB));
curPbrt.addMaterial(newMaterial);

% add material file
curPbrt.addMaterial(fullfile(s3dRootPath, 'data', 'materials', 'simpleTarget-mat.pbrt'));

%% Scene Geometry

% backdrop Depth
backDropDepth = -160;  

% calculate sphere offsets
scaleFactor = 1;

% remove default geometry
curPbrt.removeGeometry();

% add a backdrop
backDropTransform = ...
    [50 0 0 0;
    0 50 0 0 ;
    0 0 1 0;
    0 0 backDropDepth  1];
backDrop = pbrtGeometryObject('backdrop', 'Material', [], [], backDropTransform);
curPbrt.addGeometry(backDrop);

translateTransform = [1 0 0 0;
        0 1 0 0 ;
        0 0 1 0;
        0 0 0  1]; %8.87306690216     %x direction is to the right, y is into the screen, z is up
    newGeometry = pbrtGeometryObject(['sphere' ], 'grayMat', pbrtShapeObject('sphere', 'radius', 100), [], translateTransform);
    curPbrt.addGeometry(newGeometry);

sceneName = ['frontFlashInsideSphere'];
frontOi = s3dRenderOI( curPbrt, [], sceneName);
% frontFlashScene = s3dRenderScene( curPbrt, 'simpleScene');

%% Render Back flash scene

% light sources
lightBack = pbrtLightSpotObject('lightBack', [], [], [], [0 0 10], [0 0 9]);
curPbrt.removeLight();
curPbrt.addLightSource(lightBack);


% render scene
sceneName = ['backFlashInsideSphere' '.pbrt'];
% backFlashScene = s3dRenderScene( curPbrt, 'simpleScene');
backOi = s3dRenderOI( curPbrt, [], sceneName);

%% Render Depth map

% change the sampler to stratified for non-noisy depth map
samplerProp = pbrtPropertyObject();
curPbrt.sampler.setType('stratified');
curPbrt.sampler.removeProperty();
curPbrt.sampler.addProperty(pbrtPropertyObject('integer xsamples', '1'));
curPbrt.sampler.addProperty(pbrtPropertyObject('integer ysamples', '1'));
curPbrt.sampler.addProperty(pbrtPropertyObject('bool jitter', '"false"'));

% write file and render
tmpFileName = ['deleteMe'  '.pbrt'];
curPbrt.writeFile(tmpFileName);
groundTruthDepthMap = s3dRenderDepthMap(tmpFileName, 1);
figure; imagesc(groundTruthDepthMap);

frontOi = oiSet(frontOi, 'depthmap', groundTruthDepthMap);
vcAddObject(frontOi); 
backOi = oiSet(backOi, 'depthmap', groundTruthDepthMap);
vcAddObject(backOi); oiWindow;


%% front flash image processing
%load oi from file (optional)
% vcLoadObject('opticalimage', ['50mmFront.pbrt.mat']);
% oi = vcGetObject('oi');

% backOi = oi;

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

%% End