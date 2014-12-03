%% Render Scene Radiance Using pbrtObjects

clear curPbrt;
curPbrt = pbrtObject();

%camera position

%straight on
% newCamPos =    [0  0 0;
%     0   0 -1;
%     0 1.00000 0];

%angled
newCamPos =    [160  0 0;
    159   0 -1;
    0 1.00000 0];
curPbrt.camera.setPosition(newCamPos);
curPbrt.camera.lens.filmDistance = 50;
curPbrt.camera.lens.filmDiag = 70;


% curPbrt.camera.setResolution(100, 100);    %LQ mode

%uncomment to use a 2 element lens instead of a pinhole
% curPbrt.camera.setLens(fullfile(s3dRootPath, 'data', 'lens', '2ElLens50mm.pbrt'));

%sampler
sampler = curPbrt.sampler.removeProperty();
sampler.value = 512;
curPbrt.sampler.addProperty(sampler);

%backdrop Depth
% backDropDepth = -100 * scaleFactor;  %backdrop distance increases with depth of spheres
backDropDepth = -160;  
foregroundDepth = -80;
foregroundDepth2 = -70;
foregroundDepth3= -90;
%calculate sphere offsets
% xValues = linspace(-6*scaleFactor, 6*scaleFactor, 5);
% yValues = linspace(-6*scaleFactor, 6*scaleFactor, 5);
% [xOffsets yOffsets] = meshgrid(xValues, yValues); 


% lightRight = pbrtLightSpotObject('rightLight', [], [], [], inFrom, inTo);
% curPbrt.removeLight();
%lightFront = pbrtLightSpotObject('lightFront', [], [], [], [0 0 80], [0 0 -79]);
%curPbrt.addLightSource(lightFront);

curPbrt.removeLight();
% lightLeft = pbrtLightSpotObject('lightLeft', [], [], [], [20 0 -140], [19 0 -141]);
% lightLeft.setAngle(180);
% lightLeft.setDeltaAngle(180);

lightLeft = pbrtLightObject('lightLeft', 'point', [], pbrtPropertyObject('from', [5 5 -150]));
curPbrt.addLightSource(lightLeft);
% lightRight = pbrtLightSpotObject('lightRight', [], [], [], [-10 0 -100], [-9 0 -100]);
% curPbrt.addLightSource(lightRight);
% lightTop= pbrtLightSpotObject('lightTop', [], [], [], [0 10 -100], [0 9 -100]);
% curPbrt.addLightSource(lightTop);
% lightBottom= pbrtLightSpotObject('lightBottom', [], [], [], [0 -10 -100], [0 -9 -100]);
% curPbrt.addLightSource(lightBottom);


%add a new material
matRGB= [1 1 1];
newMaterial = pbrtMaterialObject('grayMat', 'matte', pbrtPropertyObject('color Kd', matRGB));
curPbrt.addMaterial(newMaterial);

%add material file
curPbrt.addMaterial(fullfile(s3dRootPath, 'data', 'materials', 'simpleTarget-mat.pbrt'));

% remove default geometry
curPbrt.removeGeometry();
%add a backdrop
backDropTransform = ...
    [50 0 0 0;
    0 50 0 0 ;
    0 0 1 0;
    0 0 backDropDepth  1];
backDrop = pbrtGeometryObject('backdrop', 'Material', [], [], backDropTransform);
curPbrt.addGeometry(backDrop);


% %add a foreground target
% foregroundTransform = ...
%     [2 0 0 0;
%     0 2 0 0 ;
%     0 0 1 0;
%     0 0 foregroundDepth  1];
% frontSquare = pbrtGeometryObject('backdrop', 'grayMat', [], [], foregroundTransform);
% curPbrt.addGeometry(frontSquare);
% 
% %add a foreground target
% foregroundTransform = ...
%     [2 0 0 0;
%     0 2 0 0 ;
%     0 0 1 0;
%     5 5 foregroundDepth2  1];
% frontSquare = pbrtGeometryObject('backdrop', 'grayMat', [], [], foregroundTransform);
% curPbrt.addGeometry(frontSquare);
% % 
% 
% %add a foreground target
% foregroundTransform = ...
%     [2 0 0 0;
%     0 2 0 0 ;
%     0 0 1 0;
%     -6 -6 foregroundDepth3  1];
% frontSquare = pbrtGeometryObject('backdrop', 'grayMat', [], [], foregroundTransform);
% curPbrt.addGeometry(frontSquare);

%serialize the x and y values for easier looping


% xOffsets = xOffsets(:);
% yOffsets = yOffsets(:);
% for i = 1:length(xOffsets)
%     %add new geoemtry
%     translateTransform = [scaleFactor 0 0 0;
%         0 scaleFactor 0 0 ;
%         0 0 scaleFactor 0;
%         xOffsets(i) yOffsets(i) sphereDepths  1]; %8.87306690216     %x direction is to the right, y is into the screen, z is up
%     newGeometry = pbrtGeometryObject(['sphere' int2str(i)], 'grayMat', pbrtShapeObject('sphere', 'radius', 1), [], translateTransform);
%     curPbrt.addGeometry(newGeometry);
% end

tmpFileName = ['deleteMe' '.pbrt'];
curPbrt.writeFile(tmpFileName);
scene = s3dRenderScene( curPbrt, 'simpleScene');

%% Render Depth map

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

scene = sceneSet(scene, 'depthmap', groundTruthDepthMap);
vcAddObject(scene); sceneWindow;