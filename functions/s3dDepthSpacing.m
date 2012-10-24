function s3dDepthSpacing(fName, inFocusDepth)
%Renders an ISET scene that has depth information.  ***NO LONGER SUPPORTED
% This probably shouldn't be a function.  There should be a tutorial script
% that calls the key functions and sets up the OI and so forth.
%
% s3d_DepthSpacing(sceneFileName, inFocusDepth)
%   sceneFileName is the .mat file containing the ISET scene.  inFocusDepth is
%   the depth in which to render in focus.  
%
% Examples:  s3d_DepthSpacing('.\column_table_goblet\ISETScene.mat',1.2)
%            s3d_DepthSpacing('.\column_table_goblet\ISETScene.mat',2.5)
%            s3d_DepthSpacing('.\dragon\ISETScene.mat', 2.1)
%            s3d_DepthSpacing('.\dragon\ISETScene.mat', 4)
%
% Renders depth image.  This is still a preliminary function.  This version does not include 
%   partial occlusion handling. 
%   The directory containing this function must be in the Matlab path for proper operation.
%
% TODO:
%   Add in better edge handling algorithms for more accurate physical
%   rendering that accounts for partial occlusion (AL).  
%   Add in support for multiple CPU's.

%% This shouldn't be in a function.  
s_initISET

%% load the scene
load(fName);
scene = sceneSet(scene,'fov',3);

%%  Make optics with a little bigger pupil

% This should be passed into the function, not created here.
oi = oiCreate;
optics = oiGet(oi,'optics');
optics = opticsSet(optics,'offaxis','cos4th');
optics = opticsSet(optics, 'otfmethod', 'custom');
optics = opticsSet(optics, 'model', 'ShiftInvariant');

fNumber = 4;
optics = opticsSet(optics,'fnumber',fNumber);

% Set the focal length large so the pupil will be large.
pupilFactor = 3;  % When set to 1, this becomes diffraction limited.
f = opticsGet(optics,'focal length');
optics = opticsSet(optics,'focal length',pupilFactor*f);

% Attach the optics to the oi and move on.
oi = oiSet(oi,'optics',optics);

%% These are the object distances for different defocus levels.

% We track depth edges over this defocus range
defocus = linspace(-1.2,0,7);

% Find the depth edges so that the range of defocus is as above, though it
% is centered around a depth of inFocusDepth.  To achieve this the
% imageDist will not be in the focal plane.
[depthEdges, imageDist, oDefocus] = oiDepthEdges(oi,defocus,inFocusDepth);

oMap  = sceneGet(scene,'depth map');
sceneDepthRange = [depthEdges(1),10];
oMap  = ieScale(oMap,sceneDepthRange(1),sceneDepthRange(2));
blurSize = 2;
supportSize = [5 5];
g = fspecial('gaussian',supportSize,blurSize); oMap = conv2(oMap,g,'same'); 
% imagesc(oMap)
scene = sceneSet(scene,'depth map',oMap);
% vcAddAndSelectObject(scene); sceneWindow;

cAberration = [];
displayFlag = 0;
[oiD, D] = oiDepthCompute(oi,scene,imageDist,depthEdges,cAberration,displayFlag);

%% Combine them and show them.
oi = oiDepthCombine(oiD,scene,depthEdges);
pupil = opticsGet(optics,'pupil radius','mm');
oi = oiSet(oi,'name',sprintf('Focus-%.1fm',inFocusDepth));
vcAddAndSelectObject(oi);  oiWindow

%%
vcSaveObject(oi);

