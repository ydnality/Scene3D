% This script renders the basic setup of 2 side-by-side light sources and a
% simple surface in between them.  There will be 2 captures.  One with each
% light source.  The hope is that we can derive some information on the
% normal vectors using this setup.  These captures are rendered using
% pbrtObjects. 


ieInit;
%% Render FrontOi
clear curPbrt;
curPbrt = pbrtObject();

%camera position
% newCamPos =    [0  0 80.0000;
%     0   0 79.0000;
%     0 1.00000 0];

newCamPos =    [0  100 80.0000;
    0   99.6 79.0000;
    0 1.00000 0];
curPbrt.camera.setPosition(newCamPos);

filmDistance = 44;
filmDiag = 50.9117;
curPbrt.camera.setLens(pbrtLensPinholeObject(filmDistance, filmDiag));  %TODO: may want to switch to a real lens later

%depths
sphereDepths = -170;  %negative is into the screen

%flash separation
flashSeparation = 6;

%scaleFactor
scaleFactor = (-sphereDepths + 80)/(80);

%backdrop Depth
% backDropDepth = -100 * scaleFactor;  %backdrop distance increases with depth of spheres
backDropDepth = sphereDepths;  %backdrop distance increases with depth of spheres

%calculate sphere offsets
xValues = linspace(-6*scaleFactor, 6*scaleFactor, 5);
yValues = linspace(-6*scaleFactor, 6*scaleFactor, 5);
[xOffsets yOffsets] = meshgrid(xValues, yValues); 


%light sources
%make a circle of lights
%use 16 different radially symmetric lights for now
lightLocationFrom = [0  100 80.0000];
lightLocationTo = [0   99.6 79.0000];
lightLocationVector = lightLocationTo - lightLocationFrom;
lightLocationVector = normvec(lightLocationVector, 'p', 2, 'dim', 2);
lightLocationBackFrom = lightLocationFrom - lightLocationVector * flashSeparation;
lightLocationBackTo = lightLocationFrom - lightLocationVector * (flashSeparation-1);
lightLocation = 80;

%reenable light cluster later.
lightCluster = false;
numLights = 16;
angleSample = linspace(0, 2*pi, numLights+1);  %+1 for the double count for 0 degrees
angleSample = angleSample(1:end-1);  %so we don't double count 0 degrees
lightClusterRadius = 2; 
lightXCoord = cos(angleSample) * lightClusterRadius;
lightYCoord = sin(angleSample) * lightClusterRadius;

% lightRight = pbrtLightSpotObject('rightLight', [], [], [], inFrom, inTo);
curPbrt.removeLight();
% curPbrt.addLightSource(lightBack);
if (lightCluster)
    for i = 1:length(lightXCoord)
        lightFront = pbrtLightSpotObject('lightFront', [], [], [], [lightXCoord(i) lightYCoord(i) lightLocation], [lightXCoord(i) lightYCoord(i)  (lightLocation - 1)]);
        curPbrt.addLightSource(lightFront);
    end
else
    lightFront = pbrtLightSpotObject('lightFront', [], [], [], lightLocationFrom, lightLocationTo);
    curPbrt.addLightSource(lightFront);
end

%add a new material
matRGB= [1 1 1];
newMaterial = pbrtMaterialObject('grayMat', 'matte', pbrtPropertyObject('color Kd', matRGB));
curPbrt.addMaterial(newMaterial);

% remove default geometry
curPbrt.removeGeometry();

%add a backdrop
backDropTransform = ...
    [100*scaleFactor 0 0 0;
    0 100*scaleFactor 0 0 ;
    0 0 1 0;
    0 0 backDropDepth  1];
backDrop = pbrtGeometryObject('backdrop', 'grayMat', [], [], backDropTransform);
curPbrt.addGeometry(backDrop);
%add floor
floorTransform = ...
    [100*scaleFactor 0 0 0;
    0 0 100*scaleFactor 0 ;
    0 -100 * scaleFactor 0 0;
    2 0 0  1];
floorPanel = pbrtGeometryObject('backdrop', 'grayMat', [], [], floorTransform);
curPbrt.addGeometry(floorPanel);

% add cone
%serialize the x and y values for easier looping
xOffsets = xOffsets(:);
yOffsets = yOffsets(:);
coneScaleFactor = 8;
translateTransform = [coneScaleFactor 0 0 0;
    0 0 * coneScaleFactor -1 * coneScaleFactor 0 ;
    0 1 * coneScaleFactor 0 *coneScaleFactor 0;
    -9.5* 2 0 -80  1]; %8.87306690216     %x direction is to the right, y is up, z is towards camera
tempShape = pbrtShapeObject('cone', 'radius', 1 * 3);
tempShape.addParameter('height', 3.5 * 3);
newGeometry = pbrtGeometryObject(['cone' int2str(i)], 'grayMat', tempShape, [], translateTransform);
curPbrt.addGeometry(newGeometry);

% add sphere
translateTransform = [1 0 0 0;
    0 1 0 0 ;
    0 0 1 0;
    18.75 * 4 6 * 4 -140  1]; %8.87306690216     %x direction is to the right, y is up, z is towards camera
tempShape = pbrtShapeObject('sphere', 'radius', 9 * 4);
newGeometry = pbrtGeometryObject(['sphere' int2str(i)], 'grayMat', tempShape, [], translateTransform);
curPbrt.addGeometry(newGeometry);

% add cylinder
  translateTransform = [1 0 0 0;
    0 0  -1  0 ;
    0 1  0  0;
    0 6 * 4 -160  1]; %8.87306690216     %x direction is to the right, y is up, z is towards camera
tempShape = pbrtShapeObject('cylinder', 'radius', 7 * 4);
tempShape.addParameter('zmin', -6 * 4);
tempShape.addParameter('zmax', 14 * 4);
newGeometry = pbrtGeometryObject(['sphere' int2str(i)], 'grayMat', tempShape, [], translateTransform);
curPbrt.addGeometry(newGeometry);

% add paraboloid

%add new geoemtry
translateTransform = [1 0 0 0;
    0 0 1 0 ;
    0 -1 0  0;
    -16* 4 20* 4 -140  1]; %8.87306690216     %x direction is to the right, y is up, z is towards camera
tempShape = pbrtShapeObject('paraboloid', 'radius', 8 * 4);
tempShape.addParameter('zmin', 0);
tempShape.addParameter('zmax', 20 * 4);
newGeometry = pbrtGeometryObject(['sphere' int2str(i)], 'grayMat', tempShape, [], translateTransform);
curPbrt.addGeometry(newGeometry);



%samplerProp = pbrtPropertyObject();
curPbrt.sampler.setType('stratified');
curPbrt.sampler.removeProperty();
curPbrt.sampler.addProperty(pbrtPropertyObject('integer xsamples', '32'));
curPbrt.sampler.addProperty(pbrtPropertyObject('integer ysamples', '32'));
curPbrt.sampler.addProperty(pbrtPropertyObject('bool jitter', '"false"'));

% curPbrt.sampler.removeProperty();
% nSamples = 64; %512;
% curPbrt.sampler.addProperty(pbrtPropertyObject('integer pixelsamples', nSamples));

sceneName = 'frontFlash';
frontOi = s3dRenderOIAndDepthMap(curPbrt, .050, sceneName);

%% render Back Oi

% lightRight = pbrtLightSpotObject('rightLight', [], [], [], inFrom, inTo);
if (lightCluster)
   for i = 1:length(lightXCoord)
       curPbrt.removeLight();
   end
else
    curPbrt.removeLight();
end

if (lightCluster)
    for i = 1:length(lightXCoord)
        lightBack = pbrtLightSpotObject('lightFront', [], [], [], [lightXCoord(i) lightYCoord(i) lightLocation + flashSeparation], [lightXCoord(i) lightYCoord(i) (lightLocation + flashSeparation - 1)]);
        curPbrt.addLightSource(lightBack);
    end
else
    lightBack = pbrtLightSpotObject('lightBack', [], [], [], lightLocationBackFrom, lightLocationBackTo);
    curPbrt.addLightSource(lightBack);
end

% curPbrt.sampler.removeProperty();
% nSamples = 64; %512;
% curPbrt.sampler.addProperty(pbrtPropertyObject('integer pixelsamples', nSamples));
% samplerProp = pbrtPropertyObject();
curPbrt.sampler.setType('stratified');
curPbrt.sampler.removeProperty();
curPbrt.sampler.addProperty(pbrtPropertyObject('integer xsamples', '32'));
curPbrt.sampler.addProperty(pbrtPropertyObject('integer ysamples', '32'));
curPbrt.sampler.addProperty(pbrtPropertyObject('bool jitter', '"false"'));

backName = 'backFlash';
backOi = s3dRenderOI(curPbrt, .050, backName);


%% get depth map
groundTruthDepthMap = oiGet(frontOi, 'depthMap');
figure; imagesc(groundTruthDepthMap);

%% front flash image processing

% uncomment this section if you wish to load oi from file
% vcLoadObject('opticalimage', ['50mmFront.pbrt.mat']);
% oi = vcGetObject('oi');
% backOi = oi;

% sensor processing
sensor = s3dProcessSensor(frontOi, 0, [400 400],0, 'analog');    %low noise, auto exposure
% sensor = s3dProcessSensor(oi, .0096, [], .03);     %high noise
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

% image processing
vciFlash = s3dProcessImage(sensor);
vcAddAndSelectObject(vciFlash); ipWindow;

%% back flash image processing

% uncomment this section if you wish to load oi from file
% vcLoadObject('opticalimage', ['50mmBack.pbrt.mat']);
% oi = vcGetObject('oi');

% sensor processing
frontFlashExpDur = sensorGet(sensor, 'expTime');
sensor = s3dProcessSensor(backOi, 0, [400 400],frontFlashExpDur, 'analog');    %low noise
% sensor = s3dProcessSensor(oi, .0096, [], .03);     %high noise
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

% image processing
vciFlashBack = s3dProcessImage(sensor);
vcAddAndSelectObject(vciFlashBack); ipWindow;

