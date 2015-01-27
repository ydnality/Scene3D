%% Volume of Linear Transforms (VOLT) ray-tracing for multi-element lenses
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
% -Check why we have a discontinuity at (0,0).  Something about the 0 always
% mapping to 0, or something ...
% -Lens element positions are all negative, and the final exit plan can be
% considered as z = 0
% -Ran into problem with z = -100, where psf was very elongated.  Run this
% again and debug.  
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
pSZ = -110:1:-60;


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
filmPosition = 159;
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

sceneFileName = fullfile(s3dRootPath, 'data', 'isetScenes',  'uniformWithDepth.mat');

if ~exist(sceneFileName,'file')
    error('No file %s\n',sceneFileName);
end

% sceneFileName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'simpleTarget', 'pinholeSceneFile.mat');

scene = load(sceneFileName);
scene = scene.scene;


%uniform depth debug
depthMap = sceneGet(scene, 'depthMap');
uniformDepth = ones(size(depthMap)) * 79;
scene = sceneSet(scene, 'depthMap', uniformDepth);

vcAddAndSelectObject(scene); sceneWindow;

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
%% Old blurring experimental code
% tic
% 
% profile on;
% 
% % Resize scene and all other related vars to the desired size(we want to oversample the film)
% unBlurredPhotons = sceneGet(scene, 'photons');
% unBlurredPhotons = imresize(unBlurredPhotons, [film.resolution(1) * 6, film.resolution(2) * 6]);
% resizedDepth = imresize(sceneGet(scene, 'depth map'), [film.resolution(1) * 6, film.resolution(2) * 6]);
% blurredPhotons = zeros(film.resolution(1), film.resolution(2), film.resolution(3));
% 
% % Declare some useful vars
% oiSize = size(unBlurredPhotons); 
% center = oiSize./2;
% sceneHFOV = sceneGet(scene, 'hfov');
% adjustedMiddleApertureRadius = 1;  %this is the size of the middle aperture
% 
% % set p the parallel pool
% 
% % originally 3580 sec
% % with parfor: 376 sec
% if (matlabpool('size') > 0)
%     matlabpool close;
% else
%     matlabpool open 8;
% end
% 
% parfor i = 1:oiSize(2)
%     i
%     film = filmC('position', [0 0 filmPosition ], ...
%     'size', [40 40], 'resolution', filmResolution,  ...
%     'wave', wave);
%     for j = 1:oiSize(1)
%         film.clear();
%         
%         %figure out rotation of pixel
%         %this will be the rotation from the positive y = 0 line in a
%         %counter-clockwise fashion
%         x = -(j - center(2));  %consider vectorizing for speed
%         y = i - center(1);
%         thetaRad = atan2(x,y);
%         thetaDeg = thetaRad/pi * 180;
%         
%         %figure out field height
%         fieldHeight = sqrt((x)^2 + (y)^2); %make into function
%         %this is only the field height with respect to the pixels on the
%         %sensor
%         
%         %convert this to the PSF location in 3 space... somehow... using
%         %the depth map and some geometry
%         wantedPSLocation = [0 0 0];
%         
%         %wantedPSLocation(2) = fieldHeight/2;  %works as a placeholder
%         %wantedPSLocation(2) = fieldHeight/oiSize(2) * resizedDepth(i,j)/filmDistance;
%         
%         %this gives the current angle with respect to optical axis, when
%         %using the radially symmetric field height (in radians)
%         currentAngle = fieldHeight/(oiSize(2)/2) * (sceneHFOV/2) * (pi/180);   % figure out FOV stuff... do we use scene FOV or oi FOV?
%         %currentDepth = 110; %assumed to be 103 for now for simplicity
%         currentDepth = resizedDepth(i,j);
%         %wantedPSLocation(3) = -currentDepth;   %old - not completely true
%         %wantedPSLocation(2) = tan(currentAngle) * currentDepth;
%         wantedPSLocation(2) = sin(currentAngle) * currentDepth;
%         wantedPSLocation(3) = -sqrt(currentDepth^2 - wantedPSLocation(2)^2); 
%         %wantedPsLocation(3) = -resizedDepth(i,j);
%         %wantedPSLocation = [0 15 -103]; %some testing with PS locations
%         
%         % --- Interpolate PSF for current point in scene ---
%         %first get the linear transform
%         LTObject = VoLTObject.interpolateAllWaves(wantedPSLocation);
%         
%         % calculate the lightfield from this particular point source
%         [inputLF]  = s3dLightFieldEntrance(wantedPSLocation, lens);   %traces rays to entrance only. 
%         
%         % Make an LT (linear transform) object and apply the LT on the inputLF
%         outputLFObject = LTObject.applyOnLF(inputLF, adjustedMiddleApertureRadius);
%         
%         % Apply linear rotation transform on LF
%         %thetaRad = theta/180 * pi;
%         rotationMatrix = [cos(thetaRad)    sin(thetaRad)      0           0;
%             -sin(thetaRad)   cos(thetaRad)      0           0;
%             0             0               cos(thetaRad)  sin(thetaRad)
%             0             0               -sin(thetaRad) cos(thetaRad)];
%         
%         rotationMatrixFull = repmat(rotationMatrix, [1 1 length(wave)]);
%         RotationObject = LTC('AInterp', rotationMatrixFull, 'wave', wave);
%         rotatedLFObject = RotationObject.applyOnLF(outputLFObject, adjustedMiddleApertureRadius);
%         
%         % Visualize PSF and phase space
%         oiI = rotatedLFObject.createOI(lens,film);
%         psfPhotons = oiGet(oiI, 'photons');
%         
%         %convert the oi of the PSF and multiply by spectral radiance of the
%         %scene, and add together, in order to blur the scene
%         
%         if(isnan(psfPhotons(:)))
%             warning('nan photons');   
%             %problem: when oiI is all 0's... for whatever reason, the min
%             %and max is not set correctly and results in nans.  Investigate
%             %this in the future;
%             psfPhotons = zeros(size(psfPhotons));   %this is a hack for now
%         end
%         
%         %weigh each channel of PSf according to the weight of the original
%         %Unblurred Image
%         illuminanceWeight = repmat(unBlurredPhotons(i,j, :), [size(psfPhotons, 1) size(psfPhotons, 2)]);
%         
%         %add blurred PSF to the existing sum
%         blurredPhotons = blurredPhotons + psfPhotons .* illuminanceWeight; 
%     end
% end
% 
% % Close the matlab pool
% if (matlabpool('size') > 0)
%     matlabpool close;
% end
% 
% % Create an oi and assign blurred photons to oi
% oi = oiCreate;
% oi = initDefaultSpectrum(oi);
% oi = oiSet(oi, 'wave', renderWave);
% oi = oiSet(oi, 'cphotons', blurredPhotons);
% oi = oiSet(oi,'hfov', hfov);
% vcAddObject(oi); oiWindow;
% 
% profile viewer;
% profile clear; 
% toc
