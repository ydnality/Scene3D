%% Render Scene Radiance Using pbrtObjects
clear curPbrt;
curPbrt = pbrtObject();

%camera position
newCamPos =    [0  0 0;
    0   0 -1;
    0 1.00000 0];
curPbrt.camera.setPosition(newCamPos);
curPbrt.camera.lens.filmDistance = 133.33;
curPbrt.camera.lens.filmDiag = 70;
curPbrt.camera.setResolution(25, 25);    %LQ mode

%uncomment to use a 2 element lens instead of a pinhole
% curPbrt.camera.setLens(fullfile(s3dRootPath, 'data', 'lens', '2ElLens50mm.pbrt'));

%sampler
sampler = curPbrt.sampler.removeProperty();
sampler.value = 512;
curPbrt.sampler.addProperty(sampler);

%backdrop Depth
backDropDepth = -100;  
foregroundDepth = -65;

%light source
lightFront = pbrtLightSpotObject('lightFront', [], [], [], [0 0 80], [0 0 -79]);
curPbrt.addLightSource(lightFront);

%add a new material
matRGB= [400 1 500 1 600 .5 700 1 ];
newMaterial = pbrtMaterialObject('grayMat', 'matte', pbrtPropertyObject('spectrum Kd', matRGB));
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

%add a foreground target
foregroundTransform = ...
    [4 0 0 0;
    0 4 0 0 ;
    0 0 1 0;
    0 0 foregroundDepth  1];
frontSquare = pbrtGeometryObject('backdrop', 'grayMat', [], [], foregroundTransform);
curPbrt.addGeometry(frontSquare);

% tmpFileName = ['deleteMe' '.pbrt'];
% curPbrt.writeFile(tmpFileName);

% radianceRenderPbrt = pbrtObject;
% radianceRenderPbrt.makeDeepCopy(curPbrt);
% scene = s3dRenderScene( radianceRenderPbrt, 'simpleScene', [], true);
% 
% %% Render Depth map
% 
% %change the sampler to stratified for non-noisy depth map
% samplerProp = pbrtPropertyObject();
% depthRenderPbrt = pbrtObject; depthRenderPbrt.makeDeepCopy(curPbrt);
% 
% %render depth map
% groundTruthDepthMap = s3dRenderDepthMap(depthRenderPbrt, 1, 'simpleScene', true);
% figure; imagesc(groundTruthDepthMap);
% 
% scene = sceneSet(scene, 'depthmap', groundTruthDepthMap);

scene = s3dRenderSceneAndDepthMap(curPbrt, 'simpleScene', true);
vcAddObject(scene); sceneWindow;