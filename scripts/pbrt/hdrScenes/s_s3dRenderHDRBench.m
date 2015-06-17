% Render the HDR bench scene
% 
% this scripts renders a pbrt scene straight up(using a .pbrt file), using
% a docker container 
%
% AL/BW Vistasoft, 2014

%% s3dRenderScene

dockerFlag = true;
inputPbrtFile = fullfile(dataPath, 'pbrtScenes','benchScene', 'sunsetHDR.pbrt'); 
scene = s3dRenderScene(inputPbrtFile, 'bench', [], dockerFlag); 

vcAddObject(scene);
sceneWindow;

%% s3dRenderSceneAndDepthMap

dockerFlag = true;
noScale    = true;
inputPbrtFile = fullfile(dataPath, 'pbrtScenes','benchScene', 'sunsetHDR.pbrt'); 

% Somehow, we need to get the associated depth file for this to run.
% Comments in the following routine.
scene = s3dRenderSceneAndDepthMap(inputPbrtFile, 'bench', dockerFlag, noScale); 

vcAddObject(scene);
sceneWindow;

%%
save('benchHDRScene','scene','inputPbrtFile');

%%