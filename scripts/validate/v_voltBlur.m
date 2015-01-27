%% Validate Volume of Linear Transforms (VOLT) ray-tracing for multi-element lenses
%
% Let's simplify
%
% This script shows how to produce a volume of linear transforms (VoLT).
% This scripts relies on a function (s3dVOLTCreateModel) that takes a lens
% object and a sampled volume of the object space and returns a structure
% of linear transforms.
%
% The linear transformation from the object position to the exit plane of
% the optics needs to account for the properties of the system apertures.
% THat is, only some of the rays from an object actually make it through
% the optical system. We discuss this issue at length in the different
% computational scripts related to this one.
%
% The idea is this:  Lens transforms of the ray from a point through the
% optics, using Snell's Law are generally not linear. But, they are locally
% linear. Thus we can produce a collection of linear transforms, one for
% each visual field position (volume), and the collection summarizes the
% overall non-linear transform.
%
% This script starts with a set of input field positions describing the
% point sources where we wish to perform complete ray-tracing.  Linear lens
% models (4x4 matrices) are computed for each these point source locations.
% These linear transforms convert the rays from each point to the output
% lightfield for that point. The sum of the light fields from all of the
% points is the complete image light field.
%
% Since the linear transforms vary slowly across the volume of points, we
% can linearly interpolate between these known estimated transforms to
% produce linear transforms for arbitrary locations of point sources.
%
%
% See also:
%    MORE SCRIPTS ILLUSTRATING STUFF ABOUT APERTURES, COMPUTATIONAL
%    EFFICIENCY, SO FORTH.
%
% Notes:
%
% AL Vistalab, 2014
%%
ieInit

%% new PS format: use spherical coordinates to specify point (this is more efficient)

pSDepth = 60:30:110;  %this is the same as before - except we use positive coordinates to be more intuitie
pSPhi = .1:5:16;  %phi will be the azimuth angle, where phi is the counter clockwise angle from the x axis
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

% Change the index of refraction for each refractive surfaces
lst     = lens.get('refractive surfaces');   % This list excludes the diaphgram
nMatrix = lens.get('n refractive surfaces'); % This list also excludes the diaphgram
nWave     = size(nMatrix,1);
nSurfaces = size(nMatrix,2);

% Add a delta to the nMatrix which changes the chromatic aberration.
% We leave the last surface unchanged because the n refers to the material
% to the right of the surface.
offset = 0.05;
deltaN  = linspace(-offset,offset,length(wave))';
for ii=1:nSurfaces
    if ~isequal(nMatrix(:,ii),ones(nWave,1))
        nMatrix(:,ii) = nMatrix(:,ii) + deltaN;
    end
end

% Put the matrix back into the n values of the refractive surfaces.
lens.set('n refractive surfaces',nMatrix);
% lens.draw


%% film (sensor) properties
% position - relative to center of final lens surface
% size - 'mm'
% wavelength samples
% filmResolution = [100 100];
filmResolution = [25 25];
filmPosition = 159;
film = filmC('position', [0 0 filmPosition ], ...
    'size', [40 40], 'resolution', filmResolution,  ...
    'wave', wave);

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
VoLTObject = VoLTC('lens', lens, 'film', film, 'fieldPositions', pSPhi, 'depths', pSDepth, 'wave', wave); 
VoLTObject.calculateMatrices();

%% Read the Scene
% If modded pbrt is NOT installed on this system, run this command to 
% load a scene file

% sceneFileName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'indestructibleObject', 'pinholeSceneFile.mat');
sceneFileName = fullfile(s3dRootPath, 'data', 'isetScenes',  'uniformWithDepth.mat');
if ~exist(sceneFileName,'file'), error('No file %s\n',sceneFileName); end

% Converted by s_convertSceneOI in ISET
scene = load(sceneFileName);
scene = scene.scene;

vcAddObject(scene); sceneWindow;

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


% Set the fov of scene so that the cropping approximately fits fov of oi (given the focal length of lens and size of film)
scene = sceneSet(scene, 'hfov', hfov);  
filmDistance = film.position(3);

%% Blur the scene - this is for a circularly symmetric lens
VoLTCameraObject = VoLTCameraC('film', film, ...
                               'VolTObject', VoLTObject, ... 
                               'scene', scene, ...
                               'lens', lens);

tic
VoLTCameraObject.blurScene([2 2]);
toc

%%
