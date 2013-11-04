%% TEST CASE - Walls straight on 
oi = s3dRenderScene('floorWallBottomBack/lambertian.pbrt');
vcAddAndSelectObject(oi);
oiWindow;

%% render depth map
% read and output depth map
% ** make sure the rendering file has a small initial aperture, and only 1
% sample per pixel!!!
depthMap = s3dRenderDepthMap('floorWallBottomBack/lambertian.pbrt', 11);
figure; imagesc(depthMap);