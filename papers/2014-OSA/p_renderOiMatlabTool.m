%% 2014 OSA Conference
%
% This script renders some chromatic aberration optical images, given PSFs
% generated from ISET.
%
% This uses the oiCompute, which is a forward calculation without including
% depth.
%
% AL Vistalab 2014
%%
s_initISET

%% Create a simple checkerboard test scene

% Let's work with a small checkerboard scene
pixPerCheck = 16;
nChecks = 2; 
scene = sceneCreate('checkerboard',pixPerCheck,nChecks);
wave  = sceneGet(scene,'wave');
scene = sceneSet(scene,'fov', 3);

% Replace the optical image into your ISET window
vcAddAndSelectObject(scene);
sceneWindow

%% Create a simple slantedbar test scene

% Let's work with a small checkerboard scene
scene = sceneCreate('slantedBar');
wave  = sceneGet(scene,'wave');
scene = sceneSet(scene,'fov', 3);

% Replace the optical image into your ISET window
vcAddAndSelectObject(scene);
sceneWindow

%% Create a simple radial lines test scene

% Let's work with a small checkerboard scene
scene = sceneCreate('radiallines');
wave  = sceneGet(scene,'wave');
scene = sceneSet(scene,'fov', 3);

% Replace the optical image into your ISET window
vcAddAndSelectObject(scene);
sceneWindow

%% Load PSF from file

% uncomment these psfFileName's if you wish to use them
%diffraction only PSF
% psfFileName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'PSFCenter_50mm_2m_f22_n401.mat');

%chromatic aberration PSF (diffraction disabled) 8m
% psfFileName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'PSFCenter_2ElLens_50mm_2mmAp_8mSensDist_n201.mat');

%chromatic aberration PSF (diffraction disabled) 2m
psfFileName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'PSFCenter_2ElLens_50mm_2mmAp_2mSensDist_n201.mat');

load(psfFileName);
oi = opticalimage;
vcAddObject(oi); oiWindow;

oiPhotons = oiGet(oi, 'photons');
numPixels = oiGet(oi,'size');
mPerSample = oiGet(oi, 'height')/numPixels(2);
oiPhotons = imresize(oiPhotons, [128 128]);   %resize to 128 because that is constricted
pixelSize = mPerSample * numPixels(2)/128 * 10^6;
umPerSample = [pixelSize pixelSize];                % Sample spacing in um

%crop the PSF so it's of managable size
% oiPhotons = oiPhotons(75-25: 75+25, :);

%% Apply PSF onto test scene

% Save the data and all the rest, in compact form 
% ieSaveSIDataFile(psf,wave,umPerSample,'customFile');
ieSaveSIDataFile(oiPhotons,wave,umPerSample,'customFile');

optics = siSynthetic('custom',oi,'customFile','deleteMe');
optics = opticsSet(optics,'model','shiftInvariant');
oi     = oiSet(oi,'optics',optics);
oi = oiCompute(scene,oi);
oi = oiSet(oi,'name','PSF Applied');

%crop image automatically
sizeDifference = oiGet(oi, 'row') - sceneGet(scene, 'row');
rows = oiGet(oi, 'row');
oi = oiCrop(oi, round([sizeDifference/2 sizeDifference/2  (rows-sizeDifference) (rows-sizeDifference)]));
vcAddAndSelectObject(oi);
oiWindow;

% imageMultiview('oi',1:4,1)

%% End