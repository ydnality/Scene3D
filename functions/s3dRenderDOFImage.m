function oi = s3dRenderDOFImageNew(scene, oi, inFocusDepth)
% Renders an ISET scene that has depth information, given the proper optics
% contained in the optical image (oi), and the desired in focus depth.  This 
% version does NOT include 
% partial occlusion handling. 
% 
% 
% Inputs:
%   scene: input scene that will be blurred.  Note that this scene MUST
%   contain a depth map!
%   oi: the optical image (which contains the lens properties)
%   inFocusDepth is: the depth in which the rendered image will be in
%   focus.
% 
% Return values:
%   oi: optical image of the depth-of-field blurred scene
% 
% Examples:  oi = s3dRenderDOFImage(scene, oi, 1.2)
%            (see user manual and s_s3dDepthSpacing.m for a tutorial on how
%            to use this function)
% 
% TODO:
%   Add in better edge handling algorithms for more accurate physical
%   rendering that accounts for partial occlusion (AL).  
%   Add in support for multiple CPU's.

%% These are the object distances for different defocus levels.

% We track depth edges over this defocus range
defocus = linspace(-1.2,0,7);

% Find the depth edges so that the range of defocus is as above, though it
% is centered around a depth of inFocusDepth.  To achieve this the
% imageDist will not be in the focal plane.
[depthEdges, imageDist, oDefocus] = oiDepthEdges(oi,defocus,inFocusDepth);

oMap  = sceneGet(scene,'depth map');  % error check in case if no depth map!!
sceneDepthRange = [depthEdges(1),10];
oMap  = ieScale(oMap,sceneDepthRange(1),sceneDepthRange(2));

% verify what this smoothing operation is doing
% blurSize = 2;
% supportSize = [5 5];
% g = fspecial('gaussian',supportSize,blurSize); oMap = conv2(oMap,g,'same'); 
% % imagesc(oMap)
scene = sceneSet(scene,'depth map',oMap);
% vcAddAndSelectObject(scene); sceneWindow;

cAberration = [];
displayFlag = 0;
%[oiD, D] =

%oiDepthCompute(oi,scene,imageDist,depthEdges,cAberration,displayFlag);
%Andy: don't need this anymore - oiDepthCompute was a wrapper that used
%s3dRenderDepthDefocus in order to calculate full-image blurs, used in the
%previous iteration.
depthCenters = oiCalculateDepthCenters(depthEdges);
%[oi, oiD, D] = s3dRenderDepthDefocusNew(scene, oi, imageDist, depthEdges, cAberration)
[oi, oiD, D] = s3dRenderDepthDefocusNew(scene, oi, imageDist, depthEdges, cAberration)

%% Combine them and return
%oi = oiDepthCombineNewPO(oiD,scene,depthEdges);     
%oi = oiDepthCombineNewPO(oiD,scene,depthCenters);     %don't need this here anymore, we are doing this inside s3dRenderDepthDefocusNew

% optics = oiGet(oi, 'optics'); 
% pupil = opticsGet(optics,'pupil radius','mm');
oi = oiSet(oi,'name',sprintf('Focus-%.1fm',inFocusDepth));
% vcAddAndSelectObject(oi);  oiWindow

%%
% vcSaveObject(oi); 

