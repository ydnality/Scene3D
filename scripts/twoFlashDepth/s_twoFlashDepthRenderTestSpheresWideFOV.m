%% Renders a test scene consisting of hemispheres with a plane to reduce shadows
%
% Andy L. Lin
%
% This script renders the basic setup of 2 side-by-side light sources and a
% simple surface in between them.  There will be 2 captures.  One with each
% light source.  The hope is that we can derive some information on the
% normal vectors using this setup.  These captures are rendered using
% pbrtObjects.  This particular script renders it using a reasonably wide
% field of view (~50 degrees) to simulate capture by a typical wide-angle
% camera.


ieInit;
%% Render FrontOi
clear curPbrt;
curPbrt = pbrtObject();

%sampler
curPbrt.sampler.setType('stratified');
curPbrt.sampler.removeProperty();
curPbrt.sampler.addProperty(pbrtPropertyObject('integer xsamples', '32'));
curPbrt.sampler.addProperty(pbrtPropertyObject('integer ysamples', '32'));
curPbrt.sampler.addProperty(pbrtPropertyObject('bool jitter', '"false"'));

%camera position
newCamPos =    [0  0 80.0000;
    0   0 79.0000;
    0 1.00000 0];
curPbrt.camera.setPosition(newCamPos);

filmDistance = 44% 74.3750; %~20 degree FOV
%35;  %this allows for a reassonably wide field of view ~42.5dFOV
filmDiag = 50.9117; %43.26;
pinholeLens = pbrtLensPinholeObject(filmDistance, filmDiag) ;
curPbrt.camera.setLens(pinholeLens)

%depths
sphereDepths = -170;  %negative is into the screen

%flash separation
flashSeparation = 6;

%scaleFactor
scaleFactor = (-sphereDepths + 80)/(80)  * 4 * .470 * 2;

%backdrop Depth
% backDropDepth = -100 * scaleFactor;  %backdrop distance increases with depth of spheres
backDropDepth = sphereDepths;  %backdrop distance increases with depth of spheres
%calculate sphere offsets
xValues = linspace(-6 * scaleFactor, 6 * scaleFactor, 5);
yValues = linspace(-6 * scaleFactor, 6 * scaleFactor, 5);
[xOffsets yOffsets] = meshgrid(xValues, yValues); 

%light sources

%make a circle of lights
%use 16 different radially symmetric lights for now
lightLocation = 80;
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
% if (lightCluster)
%     for i = 1:length(lightXCoord)
%         lightFront = pbrtLightSpotObject('lightFront', [], [], [], [lightXCoord(i) lightYCoord(i) lightLocation], [lightXCoord(i) lightYCoord(i)  (lightLocation - 1)]);
%         curPbrt.addLightSource(lightFront);
%     end
% else

%lightFront = pbrtLightSpotObject('lightFront', [], [], [], [0 0 lightLocation], [0 0 lightLocation-1]);
%curPbrt.addLightSource(lightFront);
lightFront2 = pbrtLightSpotObject('lightFront2', [], [], [], [0 0 lightLocation], [0 0 lightLocation-1]);
curPbrt.addLightSource(lightFront);
%lightFront3 = pbrtLightSpotObject('lightFront3', [], [], [], [0 0 lightLocation], [0 0 lightLocation-1]);
%curPbrt.addLightSource(lightFront);
% end


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
backDrop = pbrtGeometryObject(['sphere' int2str(i)], 'grayMat', [], [], backDropTransform);
curPbrt.addGeometry(backDrop);

%serialize the x and y values for easier looping
xOffsets = xOffsets(:);
yOffsets = yOffsets(:);
for i = 1:length(xOffsets)
    %add new geoemtry
    translateTransform = [scaleFactor 0 0 0;
        0 scaleFactor 0 0 ;
        0 0 scaleFactor 0;
        xOffsets(i) yOffsets(i) sphereDepths  1]; %8.87306690216     %x direction is to the right, y is into the screen, z is up
    newGeometry = pbrtGeometryObject(['sphere' int2str(i)], 'grayMat', pbrtShapeObject('sphere', 'radius', 1), [], translateTransform);
    curPbrt.addGeometry(newGeometry);
end

frontName = 'frontOi';
frontOi = s3dRenderOIAndDepthMap(curPbrt, .050, frontName);

%% render Back Oi

%new camera position
curPbrt.camera.setPosition(newCamPos);
curPbrt.camera.setLens(pinholeLens)

%sampler
curPbrt.sampler.setType('stratified');
curPbrt.sampler.removeProperty();
curPbrt.sampler.addProperty(pbrtPropertyObject('integer xsamples', '32'));
curPbrt.sampler.addProperty(pbrtPropertyObject('integer ysamples', '32'));
curPbrt.sampler.addProperty(pbrtPropertyObject('bool jitter', '"false"'));

%light sources
curPbrt.removeLight();
%curPbrt.removeLight();
%curPbrt.removeLight();
% if (lightCluster)
%     for i = 1:length(lightXCoord)
%         lightBack = pbrtLightSpotObject('lightFront', [], [], [], [lightXCoord(i) lightYCoord(i) lightLocation + flashSeparation], [lightXCoord(i) lightYCoord(i) (lightLocation + flashSeparation - 1)]);
%         curPbrt.addLightSource(lightBack);
%     end
% else
lightBack = pbrtLightSpotObject('lightBack', [], [], [], [0 0 lightLocation + flashSeparation], [0 0 lightLocation + flashSeparation - 1]);
curPbrt.addLightSource(lightBack);
% theta = pi/6;
% lightBack2 = pbrtLightSpotObject('lightBack2', [], [], [], [flashSeparation * sin(theta) 0 lightLocation + flashSeparation * cos(theta)], [flashSeparation * sin(theta) 0 lightLocation + flashSeparation * cos(theta) - 1]);
% curPbrt.addLightSource(lightBack2);
%lightBack3 = pbrtLightSpotObject('lightBack3', [], [], [], [-flashSeparation/2 0 lightLocation + flashSeparation/2 * sqrt(3)], [0 0 lightLocation]);
%curPbrt.addLightSource(lightBack3);

% end

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
backDrop = pbrtGeometryObject(['sphere' int2str(i)], 'grayMat', [], [], backDropTransform);
curPbrt.addGeometry(backDrop);

%serialize the x and y values for easier looping
xOffsets = xOffsets(:);
yOffsets = yOffsets(:);
for i = 1:length(xOffsets)
    %add new geoemtry
    translateTransform = [scaleFactor 0 0 0;
        0 scaleFactor 0 0 ;
        0 0 scaleFactor 0;
        xOffsets(i) yOffsets(i) sphereDepths  1]; %8.87306690216     %x direction is to the right, y is into the screen, z is up
    newGeometry = pbrtGeometryObject(['sphere' int2str(i)], 'grayMat', pbrtShapeObject('sphere', 'radius', 1), [], translateTransform);
    curPbrt.addGeometry(newGeometry);
end

backName = 'backOi';
%curPbrt.writeFile(tmpFileName);
backOi = s3dRenderOI(curPbrt, .050, backName);

%% render Back Oi 2

%new camera position
curPbrt.camera.setPosition(newCamPos);
curPbrt.camera.setLens(pinholeLens)

%sampler
curPbrt.sampler.setType('stratified');
curPbrt.sampler.removeProperty();
curPbrt.sampler.addProperty(pbrtPropertyObject('integer xsamples', '32'));
curPbrt.sampler.addProperty(pbrtPropertyObject('integer ysamples', '32'));
curPbrt.sampler.addProperty(pbrtPropertyObject('bool jitter', '"false"'));

%light sources
curPbrt.removeLight();

% lightBack = pbrtLightSpotObject('lightBack', [], [], [], [0 0 lightLocation + flashSeparation], [0 0 lightLocation ]);
% curPbrt.addLightSource(lightBack);
lightBack2 = pbrtLightSpotObject('lightBack2', [], [], [], [flashSeparation/2 0 lightLocation + flashSeparation/2 * sqrt(3)], [flashSeparation/2 0 lightLocation + flashSeparation/2 * sqrt(3) - 1]);
curPbrt.addLightSource(lightBack2);
%lightBack3 = pbrtLightSpotObject('lightBack3', [], [], [], [-flashSeparation/2 0 lightLocation + flashSeparation/2 * sqrt(3)], [0 0 lightLocation]);
%curPbrt.addLightSource(lightBack3);

% end

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
backDrop = pbrtGeometryObject(['sphere' int2str(i)], 'grayMat', [], [], backDropTransform);
curPbrt.addGeometry(backDrop);

%serialize the x and y values for easier looping
xOffsets = xOffsets(:);
yOffsets = yOffsets(:);
for i = 1:length(xOffsets)
    %add new geoemtry
    translateTransform = [scaleFactor 0 0 0;
        0 scaleFactor 0 0 ;
        0 0 scaleFactor 0;
        xOffsets(i) yOffsets(i) sphereDepths  1]; %8.87306690216     %x direction is to the right, y is into the screen, z is up
    newGeometry = pbrtGeometryObject(['sphere' int2str(i)], 'grayMat', pbrtShapeObject('sphere', 'radius', 1), [], translateTransform);
    curPbrt.addGeometry(newGeometry);
end

backName = 'backOi2';
%curPbrt.writeFile(tmpFileName);
backOi2 = s3dRenderOI(curPbrt, .050, backName);

%% get depth map
groundTruthDepthMap = oiGet(frontOi, 'depthMap');
figure; imagesc(groundTruthDepthMap);

%% front flash image processing

% uncomment this section if you wish to load oi from file
% vcLoadObject('opticalimage', ['50mmBack.pbrt.mat']);
% oi = vcGetObject('oi');

sensor = s3dProcessSensor(frontOi, 0, [400 400],0, 'analog');    %low noise, auto exposure
% sensor = s3dProcessSensor(oi, .0096, [], .03);     %high noise
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

% image processing
vciFlash = s3dProcessImage(sensor);
vcAddAndSelectObject(vciFlash); ipWindow;

%% back flash image processing

% sensor processing
frontFlashExpDur = sensorGet(sensor, 'expTime');
sensor = s3dProcessSensor(backOi, 0, [400 400],frontFlashExpDur, 'analog');    %low noise
% sensor = s3dProcessSensor(oi, .0096, [], .03);     %high noise
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

% image processing
vciFlashBack = s3dProcessImage(sensor);
vcAddAndSelectObject(vciFlashBack); ipWindow;

