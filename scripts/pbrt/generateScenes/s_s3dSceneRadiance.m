%% s_s3dSceneRadiance
%
% Render a PBRT object using pinhole optics, and save the result as an ISET
% scene.
%
% A 3D graphics image rendered using pinhole optics can be treated as the
% scene radiance.  These are useful for testing various simulations.
%
% The idea is that the pinhole in PBRT allows through only a unique ray in
% each direction.  If the image is Lambertian at every point, the stored
% image is a measure of the radiance in all directions.  (Needs to be
% explained more clearly).
% 
% The notes here show how to open a pbrt file and render the file with
% pinhole optics and then save them as a scene file (ISET).
%
% By saving these scenes and depth maps, the user can import them even
% though they don't have PBRT installed.
%
% To run these functions, though, you need to have PBRT installed.  We have
% this running on Ubuntu boxes at Stanford.  I should make a VM so that we
% can run this on Macs, too.
%
% AL Copyright Vistasoft team, 2014

%% Not properly tested!

sceneName = fullfile(s3dRootPath, 'twoFlashDepth', 'depthTargetDepths', '50mmFront.pbrt');
scene = s3dRenderScene(sceneName, 'depth target');
vcAddObject(scene); sceneWindow;

%% Depth map for depth target

% ** make sure the rendering file has a small initial aperture, and only 1
% sample per pixel rendering.

sceneName = fullfile(s3dRootPath, 'twoFlashDepth', 'depthTargetDepths','50mmDepthMap.pbrt');
depthMap = s3dRenderDepthMap(sceneName, 1);
vcNewGraphWin; imagesc(depthMap);
title('Depth-map For Depth Target (units mm)');
colorbar;

%% The metronome object image

sceneName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'indestructibleObject', 'mainPinhole.pbrt');
scene = s3dRenderScene(sceneName, 'indObj');
vcAddObject(scene); sceneWindow;
 
%% The depth map for the metronome

sceneName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'indestructibleObject', 'mainPinholeDepth.pbrt');
depthMap = s3dRenderDepthMap(sceneName, 1);
vcNewGraphWin; imagesc(depthMap);
title('Depth-map For Indestructible Object(units mm)');
colorbar;

scene = sceneSet(scene, 'depthmap', depthMap);

%% End