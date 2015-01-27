%% Volume of Linear Transforms (VOLT) ray-tracing for multi-element lenses
%
% This example calculates the apodization at each point, and takes that
% into account when blurring the image to avoid partial occlusion
% artifacts.  This is a proof-of-concept implementation, and is not
% optimized for speed, but can easily done so if implemented in C++, or
% Julia, or other efficient programming language.
%
% See also:
%    s_s3dVoLTBlurScene.m
%
% Notes:
%
% -Check why we have a discontinuity at (0,0).  Something about the 0 always
% mapping to 0, or something ...
% -Lens element positions are all negative, and the final exit plan can be
% considered as z = 0

% AL Vistalab, 2014
%%
s_initISET

%% Specify different point source positions in the object volume
%
% We compute a linear transform for each of these sample positions
%
% At this moment, we only change in field height, not yet depth.
%pSY = 0.01:.3:2;


% - causes radial marks for indObject
pSY = .1:.11:31;
pSY = [0 pSY];
pSZ = [-110 -103 -60];   %values must be monotonically increasing!!
%pSZ = -110:1:-60;


% - causes grid marks for indObject
pSY = .1:1:31;
pSY = [0 pSY];
%pSZ = [-110 -103 -60];   %values must be monotonically increasing!!
pSZ = -110:10:-60;


pSY = .1:1:31;
pSY = [0 pSY];
pSZ = [-110 -103  -90 -80 -70 -60];   %values must be monotonically increasing!!


%desired pSLocation for interpolation
wantedPSLocation = [0 15.8 -103];
%NOTE: there is something funky going at the 0 location...

%theta = -90;



%% new PS format: use spherical coordinates to specify point (this is more efficient)

pSDepth = 60:1:110;  %this is the same as before - except we use positive coordinates to be more intuitie
pSPhi = .1:1:16;  %phi will be the azimuth angle, where phi is the counter clockwise angle from the x axis
pSPhi = [0 pSPhi]; 

wantedPSLocation = [0 7 -103];  %in [degrees depth(mm)] format



%% Define the Lens and Film
lensFileName = fullfile(s3dRootPath,'data', 'lens', 'dgauss.50mm.dat');
% lensFileName = fullfile(s3dRootPath,'data', 'lens', '2ElLens.dat');

nSamples = 101;
apertureMiddleD = 1;   % mm
lens = lensC('apertureSample', [nSamples nSamples], ...
    'fileName', lensFileName, ...
    'apertureMiddleD', apertureMiddleD);
wave = lens.get('wave');

% don't change index of refraction for now... 
% Change the index of refraction for each refractive surfaces
% lst     = lens.get('refractive surfaces');   % This list excludes the diaphgram
% nMatrix = lens.get('n refractive surfaces'); % This list also excludes the diaphgram
% nWave     = size(nMatrix,1);
% nSurfaces = size(nMatrix,2);
% 
% % Add a delta to the nMatrix which changes the chromatic aberration.
% % We leave the last surface unchanged because the n refers to the material
% % to the right of the surface.
% offset = 0.05;
% deltaN  = linspace(-offset,offset,length(wave))';
% for ii=1:nSurfaces
%     if ~isequal(nMatrix(:,ii),ones(nWave,1))
%         nMatrix(:,ii) = nMatrix(:,ii) + deltaN;
%     end
% end
% 
% % Put the matrix back into the n values of the refractive surfaces.
% lens.set('n refractive surfaces',nMatrix);
% % lens.draw

% nVector = lens.surfaceArray(lst').n
% lens.set('n',nVector);
% numSurfaces = lens.get('nsurfaces');
% offset = .05;
% for i = 1:numSurfaces
%    curN = lens.surfaceArray(i).get('n');
%    
%    if (curN(1)~=0 && curN(1)~=1)
%        nVector = linspace(curN(1) - offset, curN(1) + offset, length(wave));
%        lens.surfaceArray(i).set('n', nVector);
%    end
% end

%% film (sensor) properties
% position - relative to center of final lens surface
% size - 'mm'
% wavelength samples
% filmResolution = [100 100];
filmResolution = [25 25];
%filmPosition = 159;
filmPosition = 242;
film = filmC('position', [0 0 filmPosition ], ...
    'size', [40 40], 'resolution', filmResolution,  ...
    'wave', wave);

%% Compute Snell's Law PSF at point source field height
% profile on
% 
% % For this moment, we only run at one depth
% % We set the wanted field height to the second coordinate,
% % The first coordinate is always 0.
% pointSource = wantedPSLocation;
% 
% % Compute the plenoptic pointspread that has the light field information
% % using Snell's Law.
% LF  = s3dLightField(pointSource, lens);
% oiG = LF.createOI(lens,film);
% oiG = oiSet(oiG,'name','Snell''s Law');
% vcAddAndSelectObject(oiG); oiWindow;     
% 
% % uG = plotOI(oiG,'illuminance hline',[1 135]);
% % title(sprintf(oiGet(oiG,'name')));
% profile viewer
% profile clear

%% Compute VOLT model
%
% The VOLT model consists of a collection of 4x4 A matrices that
% transforms the positions and directions of an input light-field into an
% output light-field.  Each field position has a 4x4 transform matrix. 
%
% Three collections of A matrices are calculated.  These are stored in 4x4xn
% matrices, where n is the number of input field positions.
%
%       - AComplete: the collection of complete linear transform from the
%       front-most lens element to the back-most lens element.
%       - A1stComplete: the collection of linear transforms from the
%       front-most lens element to the middle aperture. 
%       - A2ndComplete: the collection of linear transforms from the middle
%       aperture to the back-most lens element.
%
% Currently the VOLT model accomodates changes in field position ONLY.
% Depth variation will come later.
profile on
VoLTObject = VoLTC('lens', lens, 'film', film, 'fieldPositions', pSPhi, 'depths', pSDepth, 'wave', wave); 
VoLTObject.calculateMatrices();
profile viewer
profile clear


%% Read the Scene
% If modded pbrt is NOT installed on this system, run this command to 
% load a scene file

%sceneFileName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'indestructibleObject', 'pinholeSceneFile.mat');

%sceneFileName = fullfile(s3dRootPath, 'data', 'isetScenes',  'uniformWithDepth.mat');

sceneFileName = fullfile(s3dRootPath, 'data', 'pbrtScenes', 'simpleTarget', 'simpleTarget.mat');


if ~exist(sceneFileName,'file')
    error('No file %s\n',sceneFileName);
end

% sceneFileName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'simpleTarget', 'pinholeSceneFile.mat');

scene = load(sceneFileName);
scene = scene.scene;

tmpPhotons = sceneGet(scene, 'photons');
tmpPhotons = imresize(tmpPhotons, [25 25]);
scene = sceneSet(scene, 'photons', tmpPhotons);

tmpDMap = sceneGet(scene, 'depth map');
tmpDMap = imresize(tmpDMap, [25 25]);
scene = sceneSet(scene, 'depth map', tmpDMap);

% %uniform depth debug
% depthMap = sceneGet(scene, 'depthMap');
% %uniformDepth = ones(size(depthMap)) * 79;
% uniformDepth = ones([32 32]) * 100;
% uniformDepth(8:24,8:24) = 65;
% scene = sceneSet(scene, 'depthMap', uniformDepth);

vcAddAndSelectObject(scene); sceneWindow;

lens.set('aperture sample', [51 51]);  %very small amount of samples for now.


%% Some bookkeeping for scene and oi size

% Dimensions of scene data
numRows = sceneGet(scene, 'rows');
numCols = sceneGet(scene,'cols');

% Size of PSF film pixel
psfFilmPixelSize = film.size(1)/film.resolution(1);     %mm consider putting this in a function

% Film wave
renderWave = film.wave;

% Size of the scene in mm
sceneSize = sceneGet(scene, 'size');

% Use scene data as photons first
scene = sceneSet(scene,'wave', renderWave);
dM = sceneGet(scene, 'depthmap');
unBlurredPhotons = sceneGet(scene, 'photons');
vcAddAndSelectObject(scene); sceneWindow;

% Calculate approximate FOV of lens and film.  Does Michael's bbox model
%improve upon this? 
%hfov = rad2deg(2*atan2(film.size(1)/2,lens.focalLength));
% We need to take into account the magnification (wikipedia)
% alpha = 2arctan(d/2F(1 + m), where m = S2/S1 ... where 1/S2 (film distance) + 1/S1 = 1/f

% The FOV is used to calculate the uniform sampling spacing in object
% space.  We will sample a grid along the approximate FOV of the lens.  We
% will also set the scene FOV to this FOV so the crop of the sensor will
% approximately equal the FOV of the scene.
objectFocusPosition = 1/(1/lens.focalLength - 1/film.position(3)); %this is the position in the scene that will result in theoreitcal perfect focus in scene
hfov = rad2deg(2*atan2(film.size(1)/2,lens.focalLength * (1  + film.position(3)/objectFocusPosition)));
hfov = hfov * .82; %for padding - we have issues with it going off the oi right now...
%Andy: this scaling factor is very important: we need it to be just right
%such that the scene size scales perfectly to the oi pixel for pixel, or
%else there will be weird artifacts

%TODO: figure out this padding/scaling business... 

% Set the fov of scene so that the cropping approximately fits fov of oi (given the focal length of lens and size of film)
scene = sceneSet(scene, 'hfov', hfov);  
filmDistance = film.position(3);

%% Blur the scene - this is for a circularly symmetric lens

lens.apertureMiddleD = 1;
VoLTCameraObject = VoLTCameraC('film', film, ...
                               'VolTObject', VoLTObject, ... 
                               'scene', scene, ...
                               'lens', lens);

tic
%VoLTCameraObject.blurScene([2 2]);
VoLTCameraObject.blurSceneApodize([2 2]);
toc

