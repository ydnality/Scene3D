% this scripts renders a pbrt scene straight up(using a .pbrt file), using
% (1) a docker container, and also (2) using a pbrt executable.  In case 2,
% pbrt must be installed on the system and be on the system path.

%% Render Scene only with docker container
dockerFlag = true;
inputPbrtFile = fullfile(dataPath, 'pbrtScenes','indestructibleObject', 'default.pbrt'); 
scene = s3dRenderScene(inputPbrtFile, 'indestructibleObject', [], dockerFlag); 
vcAddObject(scene); sceneWindow;

%% Render scene and depth map with docker container. 
%* Note that there has to be a *_depth.pbrt file, which is identical to the
%scene file, except that it contains the following sampler:
%
% Sampler "stratified"
%    "integer xsamples" [1]
%    "integer ysamples" [1]
%    "bool jitter" "false"

dockerFlag = true;
inputPbrtFile = fullfile(dataPath, 'pbrtScenes','indestructibleObject', 'default.pbrt'); 
scene = s3dRenderSceneAndDepthMap(inputPbrtFile, 'indestructibleObject', dockerFlag); 
vcAddObject(scene); sceneWindow;