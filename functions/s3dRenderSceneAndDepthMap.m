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

if (isa(inputPbrt, 'pbrtObject'))
    %% render scene radiance
    radianceRenderPbrt = pbrtObject;
    radianceRenderPbrt.makeDeepCopy(inputPbrt);
    scene = s3dRenderScene( radianceRenderPbrt, sceneName, [], dockerFlag);

    %% Render Depth map
    %change the sampler to stratified for non-noisy depth map
    depthRenderPbrt = pbrtObject; depthRenderPbrt.makeDeepCopy(inputPbrt);
    groundTruthDepthMap = s3dRenderDepthMap(depthRenderPbrt, 1, 'simpleScene', true);
    scene = sceneSet(scene, 'depthmap', groundTruthDepthMap);
elseif (ischar(inputPbrt))
    scene = s3dRenderScene( inputPbrt, sceneName, [], dockerFlag);
    [directory, fileName, extension] = fileparts(inputPbrt);
    %depth map pbrt file must have a _depth appended to name
    depthPbrtFile = fullfile(directory, [fileName '_depth', extension]); 
    groundTruthDepthMap = s3dRenderDepthMap(depthPbrtFile, 1, 'simpleScene', true);
    scene = sceneSet(scene, 'depthmap', groundTruthDepthMap);
else
   error('invalid inputPbrt type.  Must be either a character array of the pbrt file, or a pbrtObject'); 
end
    
    
end