%% s_s3dHarmonic
%
% Perform experiments with a harmonic and different defocus measurements.
% This is good to interpret because we can see the harmonic blur as it
% comes out of the proper depth plane and the defocus changes.
%
% This script might be used to fix the algorithm for combining
% images at different depths.  For this rather harsh case, of a 
%
% (c) Stanford VISTA Team 2011
%% 
s_initISET


%%  Make optics with a little bigger pupil
oi = oiCreate;
optics = opticsCreate('default');
optics = opticsSet(optics,'offaxis','cos4th');
optics = opticsSet(optics, 'otfmethod', 'custom');
optics = opticsSet(optics, 'model', 'ShiftInvariant');

fNumber = 8;
optics = opticsSet(optics,'fnumber',fNumber);

% Set the focal length - longer means larger pupil radius
fLength = 0.050;   % 50 mm focal length
optics = opticsSet(optics,'focal length',fLength);   
% opticsGet(optics,'pupil radius')

% Find the depth edges that achieve different levels of defocus when the
% image is in the focal plane.
defocus = linspace(-1.2,-0.05,8);
depthEdges = oiDepthEdges(oi,defocus,fLength);

% Attach the optics to the oi and move on.
oi = oiSet(oi,'optics',optics);

%%  Harmonic test
parms.freq = 2; parms.contrast = 1; parms.ph = 0;
parms.ang = 0; parms.row = 128; parms.col = 128; parms.GaborFlag=0;
[scene,parms] = sceneCreate('harmonic',parms);
scene = sceneSet(scene,'fov',1);
oMap  = sceneGet(scene,'depth map');

%% Create a depth map that places different parts of the scene far enough
% away for perfect focus (depthEdges(end)) and close enough to be defocused
% by -1.2 diopters.
[r,c] = size(oMap);

% The map is also a harmonic
f = 1; 
x = 0.3*sin(2*pi*f*(1:c)/c) + 1;
d = repmat(x,r,1);
d = ieScale(d,depthEdges(1),depthEdges(end));
scene = sceneSet(scene,'depth map',d);
% vcAddAndSelectObject(scene); sceneWindow;


%% Set the best focus to be somewhere within the scene's depth range.
%
inFocusDepth = prctile(d(:),10);  % Percentile of depth that is in focus

% This defines the image plane distance that will be best for the
% inFocusDepth distance.
[depthEdges, imageDist, oDefocus] = oiDepthEdges(oi,defocus,inFocusDepth);
% This is the depth we would like to be in focus (m)

%% Create all the depth images and put them in the window

[oiD, D] = oiDepthCompute(oi,scene,imageDist,depthEdges,[],0);
for ii=1:length(oiD), vcAddAndSelectObject(oiD{ii}); end; 
oiWindow

%% The combination of depths is a problem.

% There are too many distortions with out current method.
oi = oiDepthCombine(oiD,scene,depthEdges);
vcAddAndSelectObject(oi);  oiWindow

%% End