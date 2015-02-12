%% Runs PBRT and imports it in ISET for the bench scene. 

s_initISET
%% specify the scene using pure .pbrt text file
sceneName = 'indObjAmbientOnly';

%previous input files
%('indestructibleObject/simpleReflectance-downAmbient.pbrt', .050);
%('indestructibleObject/graycard-down.pbrt', .050);
% ('desk/paintWhiteFlash.pbrt', .050);
%('desk/graycardWhiteFlash.pbrt', .050);
%('indestructibleObject/paintReflectance-downWhiteFlash.pbrt', .050);

%inputFile = fullfile(dataPath,'pbrtScenes', 'indestructibleObject/simpleReflectance-downAmbientOnly.pbrt');
%inputFile = fullfile(dataPath,'pbrtScenes', 'indestructibleObject/default.pbrt');

%% specify the scene Using pbrtObjects
clear curPbrt;
curPbrt = pbrtObject();

v1 = [-56.914787 -105.385544 13.014802];
v2 = [-56.487434 -104.481461 13.014835;];
camDir = v2 - v1;
camDir = camDir./norm(camDir, 2);
camPos = [v1 + camDir * 80 ;
    v2 + camDir* 80 ;
    -0.000013 -0.000031 1.000000;];

curPbrt.camera.setPosition(camPos);
curPbrt.camera.setLens(pbrtLensRealisticObject());

%surface integrator
curPbrt.surfaceIntegrator.setMaxDepth(1); %1 reflection

%specify camera properties
curPbrt.camera.lens.filmDistance = 15;  %90; %70; %133.33;
curPbrt.camera.lens.filmDiag = 43.75; %; 70;
%curPbrt.camera.lens.specFile = 'dgauss.50mmSA2.dat';
curPbrt.camera.lens.specFile = '2ElLens.dat';
curPbrt.camera.lens.specFile = '2ElLens13.5mm.dat';
%curPbrt.camera.lens.specFile = 'testLens.dat';


curPbrt.camera.lens.apertureDiameter = 7; % in mm
%curPbrt.camera.lens.curveRadius = 0;  %experimental
curPbrt.camera.lens.curveRadius = -13.5;  %experimental

%curPbrt.camera.setResolution(300, 200);
curPbrt.camera.setResolution(450, 300);
%curPbrt.camera.setResolution(300, 300);
%uncomment to use a 2 element lens instead of a pinhole
% curPbrt.camera.setLens(fullfile(s3dRootPath, 'data', 'lens', '2ElLens50mm.pbrt'));

% Sampler
sampler = curPbrt.sampler.removeProperty();
sampler.value = 64;
curPbrt.sampler.addProperty(sampler);

% Light sources

%spot light
spotLight = pbrtLightSpotObject();
spotLight.setName('spot');
spotLight.setSpectrum(pbrtSpectrumObject('rgb I', [1000 1000 1000]));
spotLight.setAngle(180);
spotLight.setDeltaAngle(180);
spotLight.setFrom([-142.3855 -286.2024  13.0082]);
spotLight.setTo([ -141.9582 -285.2984   13.0082]);
curPbrt.addLightSource(spotLight);

%infinite light (for diffuse lighting)
infiniteLight = pbrtLightInfiniteObject();
curPbrt.addLightSource(infiniteLight);

% Add material file
curPbrt.addMaterial(fullfile(s3dRootPath, 'data', 'pbrtScenes', 'indestructibleObject', 'default-mat.pbrt'));

% Remove default geometry
curPbrt.removeGeometry();

% Add geometry
curPbrt.addGeometry(fullfile(s3dRootPath, 'data', 'pbrtScenes', 'indestructibleObject','default-geom.pbrt'));

%% Render the oi

%oi = s3dRenderOI( inputFile, .050, sceneName);
dockerFlag = false;
oi = s3dRenderOIAndDepthMap( curPbrt, .050, sceneName, dockerFlag);


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

% %%  generate and read and output depth map
% % ** make sure the rendering file has a small initial aperture, and only 1
% % sample per pixel!
% depthMap = s3dRenderDepthMap(fullfile(dataPath, 'pbrtScenes', 'indestructibleObject/simpleReflectance-downAmbientOnly.pbrt'));
% %scene = s3dRenderScene(fullfile(dataPath, 'pbrtScenes', 'indestructibleObject', 'simpleReflectance-downAmbientOnly.pbrt'), 11);
% figure; imagesc(depthMap);