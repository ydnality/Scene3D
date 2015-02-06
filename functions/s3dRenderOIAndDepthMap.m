%This function renders an oi AND the depth map, given a pbrt object.
%Note that this function only works if the camera.lens of inputPbrt is of type
%lensPinhole. 
function oi = s3dRenderOIAndDepthMap(inputPbrt, focalLength, sceneName, dockerFlag)

if (ieNotDefined('dockerFlag'))
    dockerFlag = false;
end
if (ieNotDefined('sceneName'))
    sceneName = 'unamedScene';
end

if (ieNotDefined('focalLength'))
    focalLength = 50; %default focal length
end

if (isa(inputPbrt, 'pbrtObject'))
    %% render scene radiance
    radianceRenderPbrt = pbrtObject;
    radianceRenderPbrt.makeDeepCopy(inputPbrt);
    oi = s3dRenderOI( radianceRenderPbrt, focalLength, sceneName, dockerFlag);

    %% Render Depth map
    %change the sampler to stratified for non-noisy depth map
    depthRenderPbrt = pbrtObject; depthRenderPbrt.makeDeepCopy(inputPbrt);
    groundTruthDepthMap = s3dRenderDepthMap(depthRenderPbrt, 1, sceneName, dockerFlag);
    oi = sceneSet(oi, 'depthmap', groundTruthDepthMap);
elseif (ischar(inputPbrt))
    oi = s3dRenderOI( inputPbrt, focalLength, sceneName, dockerFlag);
    [directory, fileName, extension] = fileparts(inputPbrt);
    %depth map pbrt file must have a _depth appended to name
    depthPbrtFile = fullfile(directory, [fileName '_depth', extension]); 
    groundTruthDepthMap = s3dRenderDepthMap(depthPbrtFile, 1, sceneName, dockerFlag);
    oi = sceneSet(oi, 'depthmap', groundTruthDepthMap);
else
   error('invalid inputPbrt type.  Must be either a character array of the pbrt file, or a pbrtObject'); 
end
    
    
end