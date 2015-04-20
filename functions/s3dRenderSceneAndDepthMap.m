function scene = s3dRenderSceneAndDepthMap(inputPbrt, sceneName, dockerFlag, noScale)
% Render the scene with the depth map
%
%   scene = s3dRenderSceneAndDepthMap(inputPbrt, sceneName, dockerFlag, noScale)
%
%This function renders a scene AND the depth map, given a pbrt object.
%Note that this function only works if the camera.lens of inputPbrt is of type
%lensPinhole. 
%
if (ieNotDefined('dockerFlag')), dockerFlag = false; end
if (ieNotDefined('sceneName')),  sceneName = 'unamedScene'; end
if (ieNotDefined('noScale')),    noScale = true; end

if (isa(inputPbrt, 'pbrtObject'))
    %% render scene radiance
    radianceRenderPbrt = pbrtObject;
    radianceRenderPbrt.makeDeepCopy(inputPbrt);
    scene = s3dRenderScene( radianceRenderPbrt, sceneName, noScale, dockerFlag);

    %% Render Depth map
    %change the sampler to stratified for non-noisy depth map
    depthRenderPbrt = pbrtObject; depthRenderPbrt.makeDeepCopy(inputPbrt);
    groundTruthDepthMap = s3dRenderDepthMap(depthRenderPbrt, 1, sceneName, dockerFlag);
    scene = sceneSet(scene, 'depthmap', groundTruthDepthMap);
elseif (ischar(inputPbrt))
    %  What we want here is to convert the file inputPbrt into a pbrt
    %  object.
    
    scene = s3dRenderScene( inputPbrt, sceneName, [], dockerFlag);
    [directory, fileName, extension] = fileparts(inputPbrt);
    %depth map pbrt file must have a _depth appended to name
    depthPbrtFile = fullfile(directory, [fileName '_depth', extension]);
    if ~exist('depthPbrtFile','file')
        disp('No depth file found.  Returning scene without depth map');
        return;
    else
        numRenders = 1;
        groundTruthDepthMap = s3dRenderDepthMap(depthPbrtFile, numRenders, sceneName, dockerFlag);
        scene = sceneSet(scene, 'depthmap', groundTruthDepthMap);
    end
    
else
   error('invalid inputPbrt type.  Must be either a character array of the pbrt file, or a pbrtObject'); 
end
    
    
end