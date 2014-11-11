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

%pSLocations = [zeros(length(pSY), 1) pSY' pSZ];

%desired pSLocation for interpolation
wantedPSLocation = [0 1.7 -103];


%% Define the Lens and Film

% diffractionEnabled = false;
%turning on diffraction does NOT make sense yet since we have not modeled
%the transformation of uncertainty from the middle aperture to the end of the lens

%TODO: add diffraction into this somehow
%      Create a function lens.write() that makes a file from a lens
%      object.

%initialize and read multi-element lens from file

% both of these lenses work for this example
%lensFileName = fullfile(s3dRootPath,'data', 'lens', 'dgauss.50mm.dat');
lensFileName = fullfile(s3dRootPath,'data', 'lens', '2ElLens.dat');

nSamples = 151;
apertureMiddleD = 8;   % mm
lens = lensC('apertureSample', [nSamples nSamples], ...
    'fileName', lensFileName, ...
    'apertureMiddleD', apertureMiddleD);
wave = lens.get('wave');

% film (sensor) properties
% position - relative to center of final lens surface
% size - 'mm'
% wavelength samples
film = pbrtFilmC('position', [0 0 100 ], ...
    'size', [10 10], ...
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

%% Find centroid, and estimate maximum height/width and re-render with a 
% film at centroid, with the height and width set to twice that of the
% maximum height/width of the PSF.  Then use this new rendering form to
% find the MSE Later

% Figure out center pos by calculating the centroid of illuminance image
img = oiGet(oiG,'illuminance');

% Force to unit area and flip up/down for a point spread
img = img./sum(img(:));
img = flipud(img);
% vcNewGraphWin; mesh(img);

% Calculate the weighted centroid/center-of-mass
xSample = linspace(-film.size(1)/2, film.size(1)/2, film.resolution(1));
ySample = linspace(-film.size(2)/2, film.size(2)/2, film.resolution(2));
[filmDistanceX, filmDistanceY] = meshgrid(xSample,ySample);

distanceMatrix = sqrt(filmDistanceX.^2 + filmDistanceY.^2);
centroidX = sum(sum(img .* filmDistanceX));
centroidY = sum(sum(img .* filmDistanceY));

%% Recompute Snell's Law PSF, this time with zoomed in view
figure; [n xout] = hist(img(:), 2);

%we will pick the threshold for what is in the PSF by taking the boundary
%between the 2 histogram bins.

threshold = (xout(1) + xout(2))/2;
PSFpixels = img > threshold;
numPSFpixels = sum(PSFpixels(:));
radiusPixels = sqrt(numPSFpixels);

pixelPitch = film.size(1)/film.resolution(1);
radiusMm = radiusPixels * pixelPitch;  

zoomedSize = ones(1,2) * radiusMm * 2;

zoomedFilm = pbrtFilmC('position', [centroidX centroidY 100 ], ...
    'size', zoomedSize, ...
    'wave', wave, ...
    'resolution', [100 100]);


%Compute Snell's Law PSF at point source hield height

% For this moment, we only run at one depth
% We set the wanted field height to the second coordinate,
% The first coordinate is always 0.
pointSource = wantedPSLocation;

% Compute the plenoptic pointspread that has the light field information
% using Snell's Law.
LF  = s3dLightField(pointSource, lens);
oiSnellBig = LF.createOI(lens,zoomedFilm);
oiSnellBig = oiSet(oiSnellBig,'name','Snell''s Law');
vcAddObject(oiSnellBig); oiWindow;


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
%[AComplete A1stComplete A2ndComplete] = s3dVOLTCreateModel(lens, film, pSLocations);
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
%LTObject = LTC('wave', wave, 'AInterp', AInterp, 'A1stInterp', A1stInterp, 'A2ndInterp', A2ndInterp); 
outputLFObject = LTObject.applyOnLF(inputLF, adjustedMiddleApertureRadius);

%TODO: de-couple finding the input LF from computing ground truth

% Visualize PSF and phase space
oiI = outputLFObject.createOI(lens,film);
oiI = oiSet(oiI,'name','Light Field');
vcAddObject(oiI); oiWindow;

uI = plotOI(oiI,'illuminance hline',[1 135]);
title(sprintf(oiGet(oiI,'name')));

vcNewGraphWin;
plot(uI.pos,uI.data,'r-',uG.pos,uG.data,'b--');
grid on; xlabel('position'); ylabel('Illuminance')

vcNewGraphWin; plot(uI.data(:)/max(uI.data(:)),uG.data(:)/max(uG.data(:)),'o')
grid on; identityLine;


%% plot zoomed in VoLT PSF
adjustedMiddleApertureRadius = 4;

[~,~,inputLF]  = s3dLightField(pointSource, lens);

% Make an LT (linear transform) object and apply the LT on the inputLF
%LTObject = LTC('wave', wave, 'AInterp', AInterp, 'A1stInterp', A1stInterp, 'A2ndInterp', A2ndInterp); 
outputLFObject = LTObject.applyOnLF(inputLF, adjustedMiddleApertureRadius);

%TODO: de-couple finding the input LF from computing ground truth

% Visualize PSF and phase space
oiLFBig = outputLFObject.createOI(lens,zoomedFilm);
oiLFBig = oiSet(oiLFBig,'name','Light Field');
vcAddObject(oiLFBig); oiWindow;


%% calculate error between 2 PSF's
oiLFBigPhotons = oiGet(oiLFBig, 'photons');
oiSnellBigPhotons = oiGet(oiSnellBig, 'photons');

% normalize by wavelength
oiLFBigPhotons = oiLFBigPhotons./ ...
    repmat(sum(sum(oiLFBigPhotons,1),2), ...
    size(oiLFBigPhotons, 1), size(oiLFBigPhotons, 2));

oiSnellBigPhotons = oiSnellBigPhotons./ ...
    repmat(sum(sum(oiSnellBigPhotons,1),2), ...
    size(oiSnellBigPhotons, 1), size(oiSnellBigPhotons, 2));

squaredError = (oiLFBigPhotons - oiSnellBigPhotons).^2;
MSE = mean(squaredError(:))

%% compute the MTF of both PSF's
OTFLF = psf2otf(oiLFBigPhotons, [100 100 7]);
figure; surf(abs(fftshift(OTFLF(:,:,1))));axis square; axis tight;

OTFSnell = psf2otf(oiSnellBigPhotons, [100 100 7]);
figure; surf(abs(fftshift(OTFSnell(:,:,1))));axis square; axis tight;

% calculate error of the MTF's - and note the 50% dropoff pointow to
temp = (abs(OTFLF) - abs(OTFSnell)).^2;
MSEMTF = mean(temp(:))

%% Scrap

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



%% Calculate errors
% % Scatter plot of positions
% for ii=1:4
%     ii
%     vcNewGraphWin; plot(bOrigMasked(ii,:),bEstInterp(ii,:),'o');
%     grid on;
%     
%     meanAbsError = mean(abs(bEstInterp(ii,:) - bOrigMasked(ii,:)));
%     averageAmp = mean(abs(bOrigMasked(ii,:)));
%     meanPercentErrorSplit = meanAbsError/averageAmp * 100
% end
