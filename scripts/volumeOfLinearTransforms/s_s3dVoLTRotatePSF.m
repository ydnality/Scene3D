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
pSY = 6.01:.3:8;
pSZ = [-103 -102.75];   %values must be monotonically increasing!!

%desired pSLocation for interpolation
wantedPSLocation = [0 7.8 -103];

theta = -120;

%% Define the Lens and Film

% diffractionEnabled = false;
%turning on diffraction does NOT make sense yet since we have not modeled
%the transformation of uncertainty from the middle aperture to the end of the lens

%TODO: add diffraction into this somehow
%      Create a function lens.write() that makes a file from a lens
%      object.

%initialize and read multi-element lens from file

%lensFileName = fullfile(s3dRootPath,'data', 'lens', 'dgauss.50mm.dat');
lensFileName = fullfile(s3dRootPath,'data', 'lens', '2ElLens.dat');

nSamples = 151;
apertureMiddleD = 8;   % mm
lens = lensC('apertureSample', [nSamples nSamples], ...
    'fileName', lensFileName, ...
    'apertureMiddleD', apertureMiddleD);
wave = lens.get('wave');

%have index of refraction change with wavelength, to allow for color
%separation
%for 2 element lens
%nVector = linspace(1.5, 1.8, length(wave));
%lens.surfaceArray(1).set('n', nVector);
%lens.surfaceArray(2).set('n', nVector);

%allow for color separation with all lenses in general
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
    'size', [20 20], ...
    'wave', wave);

%% Compute Snell's Law PSF at point source hield height

% For this moment, we only run at one depth
% We set the wanted field height to the second coordinate,
% The first coordinate is always 0.
pointSource = wantedPSLocation;

% Compute the plenoptic pointspread that has the light field information
% using Snell's Law.
LF  = s3dLightField(pointSource, lens);
oiG = LF.createOI(lens,film);
oiG = oiSet(oiG,'name','Snell''s Law');
vcAddObject(oiG); oiWindow;

uG = plotOI(oiG,'illuminance hline',[1 135]);
title(sprintf(oiGet(oiG,'name')));

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
%% Obtain an A given a wanted pSLocation by linear interpolation
% A different A matrix will be calculated for each wavelengths.  We will
% loop through all the wavelengths.  The wave samples assumed are the ones
% from the lens.  These should be synchronized with everything else.  

LTObject = VoLTObject.interpolateAllWaves(wantedPSLocation);

% t1 = VoLTObject.ACollection(:,:,6,1,1)
% t2 = AInterp(:,:,1)
% vcNewGraphWin;
% plot(t1(:),t2(:),'o');
% grid on; identityLine;
%% Calculate the PSF, using the 2 A matrices and the aperture in the middle

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
vcAddObject(oiI); oiWindow;

uI = plotOI(oiI,'illuminance hline',[1 135]);
title(sprintf(oiGet(oiI,'name')));

vcNewGraphWin;
plot(uI.pos,uI.data,'r-',uG.pos,uG.data,'b--');
grid on; xlabel('position'); ylabel('Illuminance')

vcNewGraphWin; plot(uI.data(:)/max(uI.data(:)),uG.data(:)/max(uG.data(:)),'o')
grid on; identityLine;