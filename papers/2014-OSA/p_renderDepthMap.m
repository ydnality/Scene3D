%%  generate and read and output depth map for depth target
% ** make sure the rendering file has a small initial aperture, and only 1

chdir(fullfile(datapath, 'twoFlashDepth', 'depthTargetDepths'));
depthMap = s3dRenderDepthMap('50mmDepthMap.pbrt', 1);
figure; imagesc(depthMap);
title('Depth-map For Depth Target (units mm)');
colorbar;

%%  generate and read and output depth map for indestructibleObject
% ** make sure the rendering file has a small initial aperture, and only 1

chdir(fullfile(datapath, 'twoFlashDepth', 'indObject', 'pbrt'));
depthMap = s3dRenderDepthMap('50mmDepthMap.pbrt', 1);
figure; imagesc(depthMap);
title('Depth-map For Indestructible Object (units mm)');
colorbar;

