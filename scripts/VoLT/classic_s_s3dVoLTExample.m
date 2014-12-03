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
% How is this useful
%
%   This approach speeds up the forward model processing compared to a lot
%   of ray tracing through multi-element lenses.  Rather than trace each
%   ray, we can set a cone of rays for each point and apply the linear
%   transformation.
%
%   This method stores for a light field linear transformation each visual
%   field depth, field height, and wavelength.  One value of this
%   representation is that the PSF can be derived from that linear
%   transformation for any film distance.
%
%   Another advantage is that in the past we had a point spread function
%   (PSF) for every position.  We then interpolated between positions
%   (e.g., deptns and field heights) to estimate the intermediate PSF.
%
%   A third advantage is that the interpolation between positions, based on
%   the interpolation of these linear transformations, is more accurate
%   than the interpolation based on the PSFs themselves.
%
%   A fourth value is that we can recalculate the PSF as we change the
%   aperture.
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
pSY = 0.01:.3:2;
%pSZ = -102 * ones(length(pSY), 1);
pSZ = [-103 -102.75];   %values must be monotonically increasing!!

%!!weird observation: when the z spacing is too far apart, the PSF becomes
%elongated!! try to figure out why!!

%pSLocations = [zeros(length(pSY), 1) pSY' pSZ];

%desired pSLocation for interpolation
wantedPSLocation = [0 1.7 -102.9];
wantedWavelength = 550;
%% Define the Lens and Film

% diffractionEnabled = false;
%turning on diffraction does NOT make sense yet since we have not modeled
%the transformation of uncertainty from the middle aperture to the end of the lens

%TODO: add diffraction into this somehow
%      Create a function lens.write() that makes a file from a lens
%      object.

%initialize and read multi-element lens from file

% lensFileName = fullfile(s3dRootPath,'data', 'lens', 'dgauss.50mm.dat');
lensFileName = fullfile(s3dRootPath,'data', 'lens', '2ElLens.dat');

nSamples = 151;
apertureMiddleD = 8;   % mm
lens = lensC('apertureSample', [nSamples nSamples], ...
    'fileName', lensFileName, ...
    'apertureMiddleD', apertureMiddleD);

% lens.draw();
%2 element lens
% lensFile = fullfile(s3dRootPath, 'data', 'lens', '2ElLens.mat');
% import = load(lensFile,'lens');
%thickLens = import.lens;
% thickLens.apertureMiddleD = 10;

% film (sensor) properties
% position - relative to center of final lens surface
% size - 'mm'
% wavelength samples
film = pbrtFilmC('position', [0 0 100 ], ...
    'size', [10 10], ...
    'wave', 400:50:700);

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

VoLTObject = VoLTC('lens', lens, 'film', film, 'fieldPositions', pSY, 'depths', pSZ, 'wave', lens.get('wave')); 
VoLTObject.calculateMatrices();

%[AComplete A1stComplete A2ndComplete] = s3dVOLTCreateModel(lens, film, pSLocations);

%% Interpret A matrices - let's see if they vary slowly as expected

% Does A agree with MP's predict calculation from the
% method he uses?

% %% compare how A coefficients change
% for i = 1:4
%     for j = 1:4
%         ASerial = AComplete(i,j, :);
%         ASerial = ASerial(:);
%         figure; plot(pSLocations, ASerial);
%     end
% end

% Make an A movie - this is to see if the A coefficients vary slowly or
% not;

%this may go in part of the VoLTC


%**not sure if this works for multi-dimentional stuff anymore... must rework

% vcNewGraphWin; colormap(hot);
% 
% AComplete = VoLTObject.get('ACollection');
% mn = min(AComplete(:));
% mx = max(AComplete(:));
% az = -37; el = 80;
% depthIndex = 1;  %look at depth 1 for now... 
% 
% % caxis([mn mx]);
% for ii=1:size(AComplete,3)
%     surf(AComplete(:,:,ii,depthIndex)); set(gca,'zlim',[mn mx/4])
%     view(az,el);
%     shading interp
%     title(sprintf('%.2f',pSY(ii)));
%     pause(0.5);
% end
% [az el] = view;

%% Obtain an A given a wanted pSLocation by linear interpolation

[ AInterp, A1stInterp, A2ndInterp ] = VoLTObject.interpolateA(wantedPSLocation, wantedWavelength);

%% Compute ground truth LF at the WANTED Point source Field Height

% For this moment, we only run at one depth
% We set the wanted field height to the second coordinate,
% The first coordinate is always 0.
pointSource = wantedPSLocation;

% Compute the plenoptic pointspread that has the light field information
[ppsf, x, b, bMiddle, xOrig, bOrig, ppsfCamera] = s3dVOLTRTOnePoint(pointSource, film, lens);

close all;

%% Now compute the PSF

% Plot phase space and visual PSF of linear interpolation model output
ppsfCamera = s3dVOLTCreatePSFFromLF(ppsfCamera, b);

oi = ppsfCamera.oiCreate;
vcAddObject(oi); oiWindow;
plotOI(oi,'illuminance mesh linear');

%% Calculate the PSF using the INTERPOLATED A Matrix 

% The PSF for the interpolated should be close to the ground truth we just
% computed
%

% N.B. We are "cheating" in a sense because we are using the ground truth
% effective aperture information to determine which rays make it through
% the lens.
%
% FURTHER:  There is a brightness difference.  In one example, AL and I saw
% a 550 Lux vs. a 700 Lux output.  Some normalization in the interpolation
% process?  

bEstInterp = AInterp * x;

% Calculate errors
% Scatter plot of positions
for ii=1:4
    vcNewGraphWin; plot(b(ii,:),bEstInterp(ii,:),'o');
    grid on;
    
    meanAbsError = mean(abs(bEstInterp(ii,:) - b(ii,:)));
    averageAmp = mean(abs(b(ii,:)));
    meanPercentError = meanAbsError/averageAmp * 100
end

% Plot phase space and visual PSF of linear interpolation model output
ppsfCamera = s3dVOLTCreatePSFFromLF(ppsfCamera, bEstInterp)
oi = ppsfCamera.oiCreate;
vcAddObject(oi); oiWindow;
plotOI(oi,'illuminance mesh linear');


%% Calculate the PSF again, but using the 2 different A matrices (1st, 2nd)

% These matrices are again interpolated
%
% Once again, this result should not be too much different from the ground
% truth PSF.  The 2 interpolated A matrix split the lens into 2 halves, the
% first one from the front-most lens element to the aperture, and the
% second one from the aperture to the back-most lens element.
%
% However, we are "cheating" in a sense because we are using the ground
% truth effective aperture information to determine which rays make it
% through the lens.

% close all;

bEstInterp = A2ndInterp * (A1stInterp * x);

% calculate errors
% Scatter plot of positions
for ii=1:4
    ii
    vcNewGraphWin; plot(b(ii,:),bEstInterp(ii,:),'o');
    grid on;
    
    meanAbsError = mean(abs(bEstInterp(ii,:) - b(ii,:)));
    averageAmp = mean(abs(b(ii,:)));
    meanPercentErrorSplit = meanAbsError/averageAmp * 100
end

% Plot phase space and visual PSF of linear interpolation model output
ppsfCamera = s3dVOLTCreatePSFFromLF(ppsfCamera, bEstInterp)
oi = ppsfCamera.oiCreate;
vcAddObject(oi); oiWindow;
plotOI(oi,'illuminance mesh linear');

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
adjustedMiddleAperture = 2;
close all;

middleXY = ppsf.aEntranceInt.XY;
withinAperture = middleXY(:,1).^2 + middleXY(:,2).^2 <= adjustedMiddleAperture.^2;%apertureMiddleD/2;
%middleAperture = diag(middleAperture);

firstHalf = A1stInterp * xOrig;
firstHalfBlock = firstHalf(:, withinAperture); %Apply aperture to rays from the first half
bEstInterp = A2ndInterp * firstHalfBlock;

bOrigCropped = bOrig(:, withinAperture); %ground truth rays

% Calculate errors
% Scatter plot of positions
for ii=1:4
    ii
    vcNewGraphWin; plot(bOrigCropped(ii,:),bEstInterp(ii,:),'o');
    grid on;
    
    meanAbsError = mean(abs(bEstInterp(ii,:) - bOrigCropped(ii,:)));
    averageAmp = mean(abs(bOrigCropped(ii,:)));
    meanPercentErrorSplit = meanAbsError/averageAmp * 100
end

% Visualize PSF and phase space
ppsfCamera = s3dVOLTCreatePSFFromLF(ppsfCamera, bEstInterp, withinAperture)
oi = ppsfCamera.oiCreate;
vcAddObject(oi); oiWindow;
plotOI(oi,'illuminance mesh linear');

%% Interpolate a collection of A matrices for different wavelengths.  
% Next, use those to produce psf's at different wavelengths on 1 oi

% we have to concatenate the b matrix

%% End 