%%Script:  s_s3dRenderDOFImage
%
% Render a depth of field optical image.
% Illustrates how to render a depth-of-field blurred image. Reads a
% scene (with a depth map), creates an ISET optical image, optics, then 
% performs the depth of field lens simulation blurring. See manual for
% technical description of this process.  
% 
% TODO:
%   Think about better edge handling algorithms for more accurate physical
%   rendering (AL).
%   Integrate demonstration with with Maya and RenderToolBox examples


%% Initialize ISET
s_initISET

%% Set script parameters

% This is the file name of the iset scene (which MUST contain a depth map)
% to use for this depth rendering script
% fName = fullfile(s3dRootPath,'isetscenes','ISET-PianoSmall.mat');
fName = fullfile(s3dRootPath,'isetscenes','ISET-Goblet.mat');
% fName = fullfile(s3dRootPath,'isetscenes','ISET-Dragons.mat');
%  fName = fullfile(s3dRootPath,'isetscenes','ISET-livingroom.mat');
 
% fName = fullfile(s3dRootPath,'scenes', 'column_table_goblet', 'ISETScene.mat');
% fName = fullfile(s3dRootPath,'scenes', 'armadillo', 'ISETScene.mat');

% This is the depth we would like to be in focus (m)

% For ISET-Goblet
inFocusDepth = 7.0;                 %edit this value to change focal point!
%  inFocusDepth = 1.4;

% For ISET-Dragons
% inFocusDepth = 2.0;
% inFocusDepth = 4.7;
% inFocusDepth = 9.0;

% For ISET-Livingroom
% inFocusDepth = 2.2;
% inFocusDepth = 7.1;

%% load the scene
load(fName);
scene = sceneSet(scene,'fov',3);

%%  Make optics with a large pupil
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

%% Call the depth of field rendering function, and render the depth of field optical image
% oi = s3dRenderDOFImage(scene, oi, inFocusDepth)
oi = s3dRenderDOFImageNew(scene, oi, inFocusDepth)

%% Show the blurred depth of field optical image
vcAddAndSelectObject(oi);  oiWindow

%%
vcSaveObject(oi);

