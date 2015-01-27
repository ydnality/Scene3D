% this scripts renders a pbrt scene straight up(using a .pbrt file), using
% (1) a docker container, and also (2) using a pbrt executable.  In case 2,
% pbrt must be installed on the system and be on the system path.

%% Render Scene only with docker container
dockerFlag = true;
inputPbrtFile = fullfile(dataPath, 'pbrtScenes','benchScene', 'sunsetHDR.pbrt'); 
scene = s3dRenderScene(inputPbrtFile, 'bench', [], dockerFlag); 
vcAddObject(scene); sceneWindow;

