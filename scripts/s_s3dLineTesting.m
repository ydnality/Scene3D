%% OBSOLETE
%
% This doesn't use the proper ISET routines.  If you want to update, see
% s_s3dHarmonic.m
%
%
% Create a series of lines and test the defocus with distance.
%
% One way to proceed is to do this rendering on the whole scene at each
% distance and then to use the depth map to combine the whole rendering.
% 
% ieMainClose;
% clx
% ISET
% ieMainW('visible','off')

%%  Make optics with a little bigger pupil
oi = oiCreate;
optics = oiGet(oi,'optics');

% opticsGet(optics,'radius')
fNumber = 4;
optics = opticsSet(optics,'fnumber',fNumber);

% Set the focal length large so the pupil will be large.
pupilFactor = 5;
f = opticsGet(optics,'focal length');
optics = opticsSet(optics,'focal length',pupilFactor*f);

% opticsGet(optics,'radius')
optics = opticsSet(optics,'offaxis','cos4th');
optics = opticsSet(optics, 'otfmethod', 'custom');
optics = opticsSet(optics, 'model', 'ShiftInvariant');

% Attach the optics the oi and move on.
oi = oiSet(oi,'optics',optics);

%% Create an image of the line at several distances

scene = sceneCreate('lined65');
scene = sceneSet(scene,'fov',1);
sz = sceneGet(scene,'size');
fLength = opticsGet(optics,'focal length');
% imageDist = fLength*1.01;
imageDist = fLength;

% The chromatic stuff may not be working because of some of the
% calculations in OTF land.  Have a look.
cAberration = humanWaveDefocus(oiGet(oi,'wave'));
% cAberration = [];

depthEdges = [];
depthFactor = [1 1.5 2 4 6 10 100];
for ii=1:length(depthFactor)
    dMap = ones(sz(1),sz(2))*depthFactor(ii);
    scene = sceneSet(scene,'depth map',dMap);
    [oi, oiD] = s3dRenderDepthDefocus(scene,oi,imageDist,depthEdges,cAberration);
    vcAddAndSelectObject(oi);  oiWindow
end

plot(depthFactor,opticsDepthDefocus(depthFactor, optics, imageDist),'-o')
xlabel('Object depth (m)');
ylabel('Defocus (diopters)')

%% Create the image in different image planes

scene = sceneCreate('lined65');
scene = sceneSet(scene,'fov',1);

sz = sceneGet(scene,'size');
sceneDepth = 100;  % Far away
dMap = ones(sz(1),sz(2))*sceneDepth;   % Far away
scene = sceneSet(scene,'depth map',dMap);

fLength = opticsGet(optics,'focal length');

cAberration = humanWaveDefocus(oiGet(oi,'wave'));
% cAberration = [];
% The chromatic stuff may not be working because of some of the
% calculations in OTF land.  Have a look.

depthEdges = [];
planeFactor = [1 1.01, 1.02, 1.03, 1.04];
imageDist = fLength*planeFactor; 
thisD = zeros(size(imageDist));
for ii=1:length(imageDist)
    thisD(ii) = opticsDepthDefocus(100, optics, imageDist(ii));
    oi = s3dRenderDepthDefocus(scene,oi,imageDist(ii),depthEdges,cAberration);
    vcAddAndSelectObject(oi);  oiWindow
end

plot(imageDist,thisD,'-o')
xlabel('Image plane depth');
ylabel('Defocus (diopters)')
str = sprintf('Fnumber %.1f',fNumber); title(str)

%% Screw around with the piano scene
fName = fullfile(s3dRootPath,'piano_shelf','ISETSceneSmall.mat');
load(fName);
scene = sceneSet(scene,'fov',3);
oMap = sceneGet(scene,'depth map');
oMap = ieScale(oMap,1,10);

sz = sceneGet(scene,'size');
fLength = opticsGet(optics,'focal length');
% imageDist = fLength*1.01;
imageDist = fLength;

% The chromatic stuff may not be working because of some of the
% calculations in OTF land.  Have a look.
% cAberration = humanWaveDefocus(oiGet(oi,'wave'));
cAberration = [];

depthEdges = [];
depthFactor = linspace(min(oMap(:)),max(oMap(:)),5);

for ii=1:length(depthFactor)
    dMap = ones(sz(1),sz(2))*depthFactor(ii);
    scene = sceneSet(scene,'depth map',dMap);
    [oi, oiD] = s3dRenderDepthDefocus(scene,oi,imageDist,depthEdges,cAberration);
    vcAddAndSelectObject(oi);  oiWindow
end
scene = sceneSet(scene,'depth map',oMap);

plot(depthFactor,opticsDepthDefocus(depthFactor, optics, imageDist),'-o')
xlabel('Object depth (m)');
ylabel('Defocus (diopters)')
