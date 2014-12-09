%% Evaluate PSF Error using Volume of Linear Transforms (VOLT) ray-tracing
%
% See also:
%    MORE SCRIPTS ILLUSTRATING STUFF ABOUT APERTURES, COMPUTATIONAL
%    EFFICIENCY, SO FORTH.
%
% Notes:
%

% -Lens element positions are all negative, and the final exit plan can be
% considered as z = 0
%
% AL Vistalab, 2014
%%
s_initISET

%% Initialize a point,lens, film set for initial testing
s_initPLF  

% Create a polar version of pts
ptsPolar     = ptsS2P(pts);

% Sample angles and sample depths
pSTheta = 0.01:.01:.05;          % in radians
pSD     = [102.75 103 103.25];   % values must be monotonically increasing!!


%% Compute Snell's Law PSF at point source hield height

% Compute the plenoptic pointspread that has the light field information
% using Snell's Law.
LF  = s3dLightField(pts, lens);
oiG = LF.createOI(lens,film);
oiG = oiSet(oiG,'name','Snell''s Law');
vcAddObject(oiG); oiWindow;

uG = plotOI(oiG,'illuminance hline',[1 135]);
title(sprintf(oiGet(oiG,'name')));

%% Find centroid of a point oi image
%
centroid = oiGet(oiG,'centroid');

% 'distance per sample' isn't yet available in oiSet
%oiG = oiSet(oiG, 'distance per sample', 'mm', film.size(1)/film.resolution(1));
centroidmm = oiSpace(oiG,[centroid.Y,centroid.X],'mm');

% % and estimate maximum height/width and re-render with a 
% % film at centroid, with the height and width set to twice that of the
% % maximum height/width of the PSF.  Then use this new rendering form to
% % find the MSE Later
% 
% % Figure out center pos by calculating the centroid of illuminance image
% img = oiGet(oiG,'illuminance');
% 
% % Force to unit area and flip up/down for a point spread
% img = img./sum(img(:));
% img = flipud(img);
% % vcNewGraphWin; mesh(img);
% 
% % Calculate the weighted centroid/center-of-mass
% xSample = linspace(-film.size(1)/2, film.size(1)/2, film.resolution(1));
% ySample = linspace(-film.size(2)/2, film.size(2)/2, film.resolution(2));
% [filmDistanceX, filmDistanceY] = meshgrid(xSample,ySample);
% 
% distanceMatrix = sqrt(filmDistanceX.^2 + filmDistanceY.^2);
% centroidX = sum(sum(img .* filmDistanceX));
% centroidY = sum(sum(img .* filmDistanceY));

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

VoLTObject = VoLTC('lens', lens, 'film', film, 'fieldPositions', pSTheta, 'depths', pSD, 'wave', wave); 
VoLTObject.calculateMatrices();
%[AComplete A1stComplete A2ndComplete] = s3dVOLTCreateModel(lens, film, pSLocations);
%% Obtain an A given a wanted pSLocation by linear interpolation
% A different A matrix will be calculated for each wavelengths.  We will
% loop through all the wavelengths.  The wave samples assumed are the ones
% from the lens.  These should be synchronized with everything else.  

LTObject = VoLTObject.interpolateAllWaves(wantedPSLocPolar);

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
film.clear();
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
zoomedFilm.clear();
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
