%This function renders a scene AND the depth map, given a pbrt object.
%Note that this function only works if the camera.lens of inputPbrt is of type
%lensPinhole. 
function scene = s3dRenderSceneAndDepthMap(inputPbrt, sceneName, dockerFlag)

if (ieNotDefined('dockerFlag'))
    dockerFlag = false;
end
if (ieNotDefined('sceneName'))
    sceneName = 'unamedScene';
end

%% render scene radiance
radianceRenderPbrt = pbrtObject;
radianceRenderPbrt.makeDeepCopy(inputPbrt);
scene = s3dRenderScene( radianceRenderPbrt, sceneName, [], dockerFlag);

%% Render Depth map
%change the sampler to stratified for non-noisy depth map
depthRenderPbrt = pbrtObject; depthRenderPbrt.makeDeepCopy(inputPbrt);

%render depth map
groundTruthDepthMap = s3dRenderDepthMap(depthRenderPbrt, 1, 'simpleScene', true);
scene = sceneSet(scene, 'depthmap', groundTruthDepthMap);
%vcAddObject(scene); sceneWindow;
end