%This function renders a scene AND the depth map, given a pbrt object.
%Note that this function only works if the camera.lens of inputPbrt is of type
%lensPinhole. 
function scene = s3dRenderSceneAndDepthMap(inputPbrt, sceneName)

%places all generated pbrt files into the generatedPbrtFiles directory
fullfname = fullfile(dataPath, 'generatedPbrtFiles', [sceneName '.pbrt']);
inputPbrt.writeFile(fullfname);

scene = s3dRenderScene(fullfname, sceneName);

%% depth map - consider putting this into it's own function

%inputs - pbrt object, output, rendered scene, and rendered depth map
%the lens must be a pinhole camera, or a lens with a pinhole aperture
inputPbrt.sampler.setType('stratified');
inputPbrt.sampler.removeProperty();
inputPbrt.sampler.addProperty(pbrtPropertyObject('bool jitter', '"false"'));
inputPbrt.sampler.addProperty(pbrtPropertyObject('integer xsamples', 1));
inputPbrt.sampler.addProperty(pbrtPropertyObject('integer ysamples', 1));

fullfnameDM = fullfile(dataPath, 'generatedPbrtFiles', [sceneName 'DM.pbrt']);
inputPbrt.writeFile(fullfnameDM);
depthMap = s3dRenderDepthMap(fullfnameDM, 1);  %this is in mm
scene = sceneSet(scene, 'depthMap', depthMap);

end