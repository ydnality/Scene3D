%% p_renderDepthMap
%
% This script generates the depth maps for 2 scenes.  The first scene is a
% test scene featuring many depth targets.  The second scene is the
% indestructible Object scene.
%
% AL

error('Obsolete.  See s_s3dSceneRadiance');

%
datapath = fullfile(s3dRootPath,'data');

%%  generate and read and output depth map for depth target
% ** make sure the rendering file has a small initial aperture, and only 1
% sample per pixel rendering.

chdir(fullfile(datapath, 'twoFlashDepth', 'depthTargetDepths'));
depthMap = s3dRenderDepthMap('50mmDepthMap.pbrt', 1);
figure; imagesc(depthMap);
title('Depth-map For Depth Target (units mm)');
colorbar;

%%  generate and read and output depth map for indestructibleObject
% ** make sure the rendering file has a small initial aperture, and only 1
% sample per pixel rendering.

chdir(fullfile(datapath, 'twoFlashDepth', 'indObject', 'pbrt'));
depthMap = s3dRenderDepthMap('50mmDepthMap.pbrt', 1);
figure; imagesc(depthMap);
title('Depth-map For Indestructible Object (units mm)');
colorbar;

