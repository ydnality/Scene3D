%% Render Scene Radiance for some simple squares using pbrtObjects
%
%  

%%
clear curPbrt;
curPbrt = pbrtObject();

%% Set camera position
newCamPos =    [0  0 0;
    0   0 -1;
    0 1.00000 0];
curPbrt.camera.setPosition(newCamPos);
curPbrt.camera.lens.filmDistance = 133.33;
curPbrt.camera.lens.filmDiag = 70;
curPbrt.camera.setResolution(50, 50);    %LQ mode

%uncomment to use a 2 element lens instead of a pinhole
% curPbrt.camera.setLens(fullfile(s3dRootPath, 'data', 'lens', '2ElLens50mm.pbrt'));

% Sampler determines how many rays are traced per pixel
sampler = curPbrt.sampler.removeProperty();  % Eliminate the default
sampler.value = 512;                         % Replace with this
curPbrt.sampler.addProperty(sampler);

% Remove default geometry from the object
% removeGeometry() means remove top of the stack.
% removeGeometry(i) means remove the ith element of the stack
% Could have a replaceGeometry(i) which would be remove/add 
curPbrt.removeGeometry();

%% Set up the background
backDropDepth = -100;  % Millimeters
% Sets up the scale and translation for the background
backDropTransform = ...
    [50 0 0 0;
    0 50 0 0 ;
    0 0 1 0;
    0 0 backDropDepth  1];
% inName, inMaterial, inTriMesh, inPoints, inTransform
% The default mesh is a simple planar object.
backDrop = pbrtGeometryObject('backdrop', 'Material', [], [], backDropTransform);
curPbrt.addGeometry(backDrop);

%% Light source
% (inName, inSpectrum, inConeAngle, inDeltaAngle, inFrom, inTo)
% Positive values are behind the lens and negative are in object space
% This one is straight on.s
curPbrt.removeLight();
inSpectrum = pbrtSpectrumObject('spectrum I', [400 0 500 .3 600 1 700 1]);
lgt = pbrtLightSpotObject('lightFront', inSpectrum, [], [], [0 0 80], [0 0 79]);
curPbrt.addLightSource(lgt);

% Second light source
inSpectrum = pbrtSpectrumObject('spectrum', [400 1 500 .1 600 .1 700 .1]);
lgt = pbrtLightSpotObject('lightUpwards', inSpectrum, [], [], [0 75 75], [0 70 70]);
curPbrt.addLightSource(lgt);

%% Add a new material
curPbrt.addMaterial(fullfile(s3dRootPath, 'data', 'materials', 'simpleTarget-mat.pbrt'));

%% A couple of foreground squares
%  We will put them at different depths and givem slightly different
%  material properties

% Add a foreground target.  Notice the x,y scales are smaller than the
% background target.
foregroundDepth = -65; xSize = 4; ySize = 4;
foregroundTransform = ...
    [xSize 0 0 0;
    0 ySize 0 0 ;
    0 0 1 0;
    0 0 foregroundDepth  1];

% Diffuse reflectance (Kd) at each wavelength
matRGB= [400 1 500 1 600 .5 700 1 ];
newMaterial = pbrtMaterialObject('cyanMat', 'matte', pbrtPropertyObject('spectrum Kd', matRGB));
curPbrt.addMaterial(newMaterial);

s1 = pbrtGeometryObject('square1', 'cyanMat', [], [], foregroundTransform);
curPbrt.addGeometry(s1);

% % Second square
foregroundDepth = -40; xSize = 1.5; ySize = 1.5;
foregroundTransform = ...
    [xSize 0 0 0;
    0 ySize 0 0 ;
    0 0 1 0;
    0 0 foregroundDepth  1];
matRGB= [400 1 500 .3 600 .3 700 0 ];
newMaterial = pbrtMaterialObject('blueMat', 'matte', pbrtPropertyObject('spectrum Kd', matRGB));
curPbrt.addMaterial(newMaterial);

s2 = pbrtGeometryObject('square2', 'blueMat', [], [], foregroundTransform);
curPbrt.addGeometry(s2);

%% Render scene and depth map

% Calls docker container to run PBRT
scene = s3dRenderSceneAndDepthMap(curPbrt, 'simpleScene', true);
vcAddObject(scene); sceneWindow;

%%

