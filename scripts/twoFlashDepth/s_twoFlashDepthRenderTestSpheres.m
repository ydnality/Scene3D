% This script renders the basic setup of 2 side-by-side light sources and a
% simple surface in between them.  There will be 2 captures.  One with each
% light source.  The hope is that we can derive some information on the
% normal vectors using this setup.  These captures are rendered using
% pbrtObjects. 




% s_initISET
%% Render FrontOi
clear curPbrt;
curPbrt = pbrtObject();

%camera position
newCamPos =    [0  0 80.0000;
    0   0 79.0000;
    0 1.00000 0];
curPbrt.camera.setPosition(newCamPos);

%depths
sphereDepths = 30;  %negative is into the screen

%flash separation
flashSeparation = 50;

%scaleFactor
scaleFactor = (-sphereDepths + 80)/(80);

%backdrop Depth
backDropDepth = -100 * scaleFactor;  %backdrop distance increases with depth of spheres

%calculate sphere offsets
xValues = linspace(-6*scaleFactor, 6*scaleFactor, 5);
yValues = linspace(-6*scaleFactor, 6*scaleFactor, 5);
[xOffsets yOffsets] = meshgrid(xValues, yValues); 

%light sources
lightFront = pbrtLightSpotObject('lightFront', [], [], [], [0 0 80], [0 0 79]);
lightBack = pbrtLightSpotObject('lightBack', [], [], [], [0 0 80 + flashSeparation], [0 0 79 + flashSeparation]);
% lightRight = pbrtLightSpotObject('rightLight', [], [], [], inFrom, inTo);
curPbrt.removeLight();
% curPbrt.addLightSource(lightBack);
curPbrt.addLightSource(lightFront);

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

tmpFileName = ['deleteMe' '.pbrt'];
curPbrt.writeFile(tmpFileName);
frontOi = s3dRenderOI(curPbrt, .050, tmpFileName);

%% render Back Oi
clear curPbrt;
curPbrt = pbrtObject();

% %camera position
% newCamPos =    [0  0 80.0000;
%     0   0 79.0000;
%     0 1.00000 0];

% 
% %depths
% sphereDepths = -200;  %negative is into the screen
% 
% %flash separation
% flashSeparation = 50;
% 
% %scaleFactor
% scaleFactor = (-sphereDepths + 80)/(80);
% 
% %backdrop Depth
% backDropDepth = -100 * scaleFactor;  %backdrop distance increases with depth of spheres

%new camera position
curPbrt.camera.setPosition(newCamPos);

%calculate sphere offsets
xValues = linspace(-6*scaleFactor, 6*scaleFactor, 5);
yValues = linspace(-6*scaleFactor, 6*scaleFactor, 5);
[xOffsets yOffsets] = meshgrid(xValues, yValues); 

%light sources
lightFront = pbrtLightSpotObject('lightFront', [], [], [], [0 0 80], [0 0 79]);
lightBack = pbrtLightSpotObject('lightBack', [], [], [], [0 0 80 + flashSeparation], [0 0 79 + flashSeparation]);
% lightRight = pbrtLightSpotObject('rightLight', [], [], [], inFrom, inTo);
curPbrt.removeLight();
curPbrt.addLightSource(lightBack);
% curPbrt.addLightSource(lightFront);

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

tmpFileName = ['deleteMe' '.pbrt'];
curPbrt.writeFile(tmpFileName);
backOi = s3dRenderOI(curPbrt, .050, tmpFileName);


%% render depth map

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

