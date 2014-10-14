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

pSY = 15:.3:17;
pSZ = [-103 -102.75];   %values must be monotonically increasing!!

%desired pSLocation for interpolation
wantedPSLocation = [0 15.8 -103];

theta = -90;


%% Define the Lens and Film

%lensFileName = fullfile(s3dRootPath,'data', 'lens', 'dgauss.50mm.dat');
lensFileName = fullfile(s3dRootPath,'data', 'lens', '2ElLens.dat');

nSamples = 151;
apertureMiddleD = 8;   % mm
lens = lensC('apertureSample', [nSamples nSamples], ...
    'fileName', lensFileName, ...
    'apertureMiddleD', apertureMiddleD);
wave = lens.get('wave');

% Change the index of refraction for the surfaces
%TODO: put this into a function for lens

numSurfaces = lens.get('nsurfaces');
offset = .05;
for i = 1:numSurfaces
   curN = lens.surfaceArray(i).get('n');
   
   if (curN(1)~=0 && curN(1)~=1)
       nVector = linspace(curN(1) - offset, curN(1) + offset, length(wave));
       lens.surfaceArray(i).set('n', nVector);
   end
end

% film (sensor) properties
% position - relative to center of final lens surface
% size - 'mm'
% wavelength samples
film = pbrtFilmC('position', [0 0 100 ], ...
    'size', [40 40], 'resolution', [50 50],  ...
    'wave', wave);

%% Compute Snell's Law PSF at point source field height

% For this moment, we only run at one depth
% We set the wanted field height to the second coordinate,
% The first coordinate is always 0.
pointSource = wantedPSLocation;

% Compute the plenoptic pointspread that has the light field information
% using Snell's Law.
LF  = s3dLightField(pointSource, lens);
oiG = LF.createOI(lens,film);
oiG = oiSet(oiG,'name','Snell''s Law');
vcAddAndSelectObject(oiG); oiWindow;     

% uG = plotOI(oiG,'illuminance hline',[1 135]);
% title(sprintf(oiGet(oiG,'name')));

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

VoLTObject = VoLTC('lens', lens, 'film', film, 'fieldPositions', pSY, 'depths', pSZ, 'wave', wave); 
VoLTObject.calculateMatrices();



%% Read the Scene
% If modded pbrt is NOT installed on this system, run this command to 
% load a scene file

sceneFileName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'indestructibleObject', 'pinholeSceneFile.mat');
% sceneFileName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'simpleTarget', 'pinholeSceneFile.mat');

scene = load(sceneFileName);
scene = scene.scene;
vcAddAndSelectObject(scene); sceneWindow;

%% Determine the size of the output image

% Dimensions of scene data
numRows = sceneGet(scene, 'rows');
numCols = sceneGet(scene,'cols');

% Size of PSF film pixel
%smallPixelSize = PSFStructure.film.size(1)/size(PSFStructure.film.image,1);     %mm
psfFilmPixelSize = film.size(1)/film.resolution(1);     %mm consider putting this in a function

% Size of scene sample pixel
%largePixelSize = 70/sqrt(2)/numRows;  %mm   
renderWave = film.wave;

%should match the PSFStructure film Z to be consistent
filmDistance = film.position(3); 

%70/sqrt(2) represents row size/diagonal  TODO: find a way to automate this
%This size should correspond to the size sensor you wish to use for the
%forward calculation.  

% Amount that PSF film needs to be scaled to be equivalent to scene sample
%size
%scaleFactor = smallPixelSize/largePixelSize;
sceneSize = sceneGet(scene, 'size');
scaleFactor = film.resolution(1)/sceneSize(1); %conversion factor from scene scale to film scale

% Create an oi
oi = oiCreate;
oi = initDefaultSpectrum(oi);

% Use scene data as photons first
scene = sceneSet(scene,'wave', renderWave);
dM = sceneGet(scene, 'depthmap');
photons = sceneGet(scene, 'photons');
oi = oiSet(oi, 'wave', renderWave);
oi = oiSet(oi, 'cphotons', photons);


%experimental  - this is copied from psfCamera.oiCreate.  Is this correct?
%What are some non ideal behavior?   do we take from size(1), or size(2)?
%for horizontal FOV?  
hfov = rad2deg(2*atan2(film.size(1)/2,lens.focalLength));
oi = oiSet(oi,'hfov', hfov);
filmDistance = film.position(3);


% Pad the oi for future processing
% padAmount = (round(size(PSFStructure.film.image,1) * scaleFactor) * 2);
% oi = oiPad(oi, [padAmount padAmount]);

% Obtain unblurred image, to be used later
unBlurredPhotons = oiGet(oi, 'photons');
photonSum = zeros(size(oiGet(oi,'photons')));

vcAddObject(oi); oiWindow;

%% Blur the scene - this is for a circularly symmetric lens

%resize scene to the size of film
unBlurredPhotons = oiGet(oi, 'photons');
unBlurredPhotons = imresize(unBlurredPhotons, [film.resolution(1), film.resolution(2)]);
resizedDepth = imresize(sceneGet(scene, 'depth map'), [film.resolution(1), film.resolution(2)]);
blurredPhotons = zeros(size(unBlurredPhotons));

%loop through each pixel
%for each pixel...
oiSize = size(unBlurredPhotons); 
center = oiSize./2;

for i = 1:oiSize(1)
    i
    for j = 1:oiSize(2)
        
        %figure out rotation of pixel
        %this will be the rotation from the positive y = 0 line in a
        %counter-clockwise fashion
        x = i - center(2);  %consider vectorizing for speed
        y = j - center(1);
        thetaRad = atan(y/x);
        thetaDeg = thetaRad/pi * 180;
        
        %figure out field height
        fieldHeight = sqrt((x)^2 + (y)^2); %make into function
        %this is only the field height with respect to the pixels on the
        %sensor
        
        %convert this to the PSF location in 3 space... somehow... using
        %the depth map and some geometry
        
        
        wantedPSLocation = [0 0 0];
        wantedPSLocation(3) = -103;
        %wantedPSLocation(2) = fieldHeight/oiSize(2) * resizedDepth(i,j)/filmDistance;
        %wantedPsLocation(3) = -resizedDepth(i,j);
        
        wantedPSLocation = [0 15 -103];
        
        
        % --- figure out PSF for current point in scene ---
        
        %first get the linear transform
        LTObject = VoLTObject.interpolateAllWaves(wantedPSLocation);
        
        %calculate the lightfield from this particular point source
        adjustedMiddleApertureRadius = 4;
        [inputLF]  = s3dLightFieldEntrance(wantedPSLocation, lens);   %traces rays to entrance only. 
        
        % Make an LT (linear transform) object and apply the LT on the inputLF
        outputLFObject = LTObject.applyOnLF(inputLF, adjustedMiddleApertureRadius);
        
        % Apply linear rotation transform on LF
        %thetaRad = theta/180 * pi;
        rotationMatrix = [cos(thetaRad)    sin(thetaRad)      0           0;
            -sin(thetaRad)   cos(thetaRad)      0           0;
            0             0               cos(thetaRad)  sin(thetaRad)
            0             0               -sin(thetaRad) cos(thetaRad)];
        
        
        rotationMatrixFull = repmat(rotationMatrix, [1 1 length(wave)]);
        RotationObject = LTC('AInterp', rotationMatrixFull, 'wave', wave);
        rotatedLFObject = RotationObject.applyOnLF(outputLFObject, adjustedMiddleApertureRadius);
        
        % Visualize PSF and phase space
        oiI = rotatedLFObject.createOI(lens,film);  %has 7 channels for some reason
        psfPhotons = oiGet(oiI, 'photons');
        
        
        %convert the oi of the PSF and multiply by spectral radiance of the
        %scene, and add together, in order to blur the scene
        
        illuminanceWeight = repmat(unBlurredPhotons(i,j, :), [size(psfPhotons, 1) size(psfPhotons, 2)]);
        blurredPhotons = blurredPhotons + psfPhotons .* illuminanceWeight; 

    end
end

%show blurredPhotons somehow
maxVal = max(max(blurredPhotons(:,:,1)));
figure; imshow(blurredPhotons(:,:,1)/maxVal);

%% 
% Obtain an A given a wanted pSLocation by linear interpolation
% A different A matrix will be calculated for each wavelengths.  We will
% loop through all the wavelengths.  The wave samples assumed are the ones
% from the lens.  These should be synchronized with everything else.  

LTObject = VoLTObject.interpolateAllWaves(wantedPSLocation);

% Calculate the PSF, using the 2 A matrices and the aperture in the middle

% This illustrates that the aperture can be changed quickly in the
% simulation, without recomputing the VOLT class
%
% Once again, this result should not be too much different from the ground
% truth PSF, unless adjustedMiddleAperture was changed.
%
% This experiment demonstrates the flexibility of the transform method.  We
% can quickly produce PSFs at arbitrary field heights, and aperture shapes
% and sizes, given the VOLT model.
%
% We are no longer "cheating" in this experiemnt because we are using the
% assumption that ONLY the middle aperture will constrict light flow in
% this system.  For middle apertures that are smaller than the other
% apertures, this assumption should be valid.

% Break this out into a separate couple of scripts that illustrate changing
% the properties of the aperture

%*** change this parameter to change the size of the middle aperture for the
%lens


% new film...
film = pbrtFilmC('position', [0 0 100 ], ...
    'size', [40 40], ...
    'wave', wave);


adjustedMiddleApertureRadius = 4;

[~,~,inputLF]  = s3dLightField(pointSource, lens);

% Make an LT (linear transform) object and apply the LT on the inputLF
outputLFObject = LTObject.applyOnLF(inputLF, adjustedMiddleApertureRadius);

% Apply linear rotation transform on LF
thetaRad = theta/180 * pi;
rotationMatrix = [cos(thetaRad)    sin(thetaRad)      0           0;
                  -sin(thetaRad)   cos(thetaRad)      0           0;
                  0             0               cos(thetaRad)  sin(thetaRad)
                  0             0               -sin(thetaRad) cos(thetaRad)];
              

rotationMatrixFull = repmat(rotationMatrix, [1 1 length(wave)]);              
RotationObject = LTC('AInterp', rotationMatrixFull, 'wave', wave);
rotatedLFObject = RotationObject.applyOnLF(outputLFObject, adjustedMiddleApertureRadius);
              
% Visualize PSF and phase space
oiI = rotatedLFObject.createOI(lens,film);
oiI = oiSet(oiI,'name','Light Field');
vcAddAndSelectObject(oiI); oiWindow;

theta = -310;
thetaRad = theta/180 * pi;
rotationMatrix = [cos(thetaRad)    sin(thetaRad)      0           0;
                  -sin(thetaRad)   cos(thetaRad)      0           0;
                  0             0               cos(thetaRad)  sin(thetaRad)
                  0             0               -sin(thetaRad) cos(thetaRad)];
              

rotationMatrixFull = repmat(rotationMatrix, [1 1 length(wave)]);              
RotationObject = LTC('AInterp', rotationMatrixFull, 'wave', wave);
rotatedLFObject = RotationObject.applyOnLF(outputLFObject, adjustedMiddleApertureRadius);
              




% Visualize PSF and phase space
oiI = rotatedLFObject.createOI(lens,film);
oiI = oiSet(oiI,'name','Light Field');
vcAddAndSelectObject(oiI); oiWindow;

% uI = plotOI(oiI,'illuminance hline',[1 135]);
% title(sprintf(oiGet(oiI,'name')));
% 
% vcNewGraphWin;
% plot(uI.pos,uI.data,'r-',uG.pos,uG.data,'b--');
% grid on; xlabel('position'); ylabel('Illuminance')

% vcNewGraphWin; plot(uI.data(:)/max(uI.data(:)),uG.data(:)/max(uG.data(:)),'o')
% grid on; identityLine;
