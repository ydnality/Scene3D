%% Volume of Linear Transforms (VOLT) ray-tracing for multi-element lenses 
%
% This script produces a volume of linear transforms (VOLT) to help speed
% up ray-tracing.  Although lens transforms using Snell's Law are generally
% not linear, it can be shown that they are approximately locally linear.
% Thus we can produce a collection of linear transforms to summarize the
% generally non-linear transform.
%
% Input field positions are given for the point sources where we wish to
% perform complete ray-tracing.  Linear lens models (4x4 matrices) are
% computed for these point source locations that transform the input light
% field to the output lightfield.
%
% Since the linear transforms were shown to vary slowly, we can linearly
% interpolate between these known estimated transforms to produce linear
% transforms for arbitrary locations of point sources.  
%
% Notes:
%
% -Check why we have a discontinuity at (0,0).  Something about the 0 always
% mapping to 0, or something ...
% -Lens element positions are all negative, and the final exit plan can be
% considered as z = 0
%
% AL Vistalab, 2014
%%
s_initISET

%% Specify different point source positions to compute linear transform
pSY = 0.01:.3:2;
pSZ = -102 * ones(length(pSY), 1);

pSLocations = [zeros(length(pSY), 1) pSY' pSZ];

%desired pSLocation for interpolation
wantedPSLocation = 1.7;

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
% output light-field.  Each field position has a 4x4 transform matrix. 3
% collections of A matrices are calculated.  These are stored in 4x4xn
% matrices, where n is the number of input field positions.
%       -AComplete: the collection of complete linear transform from the
%       front-most lens element to the back-most lens element.
%       -A1stComplete: the collection of linear transforms from the
%       front-most lens element to the middle aperture. -A2ndComplete: the
%       collection of linear transforms from the middle aperture to the
%       back-most lens element.
%
% Currently the VOLT model accomodates changes in field position ONLY.
% Depth variation will come later.

[AComplete A1stComplete A2ndComplete] = s3dVOLTCreateModel(lens, film, pSLocations);

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
vcNewGraphWin; colormap(hot);
mn = min(AComplete(:));
mx = max(AComplete(:));
az = -37; el = 80;
% caxis([mn mx]);
for ii=1:size(AComplete,3)
    surf(AComplete(:,:,ii)); set(gca,'zlim',[mn mx/4])
    view(az,el);
    shading interp
    title(sprintf('%.2f',pSY(ii)));
    pause(0.5);
end
% [az el] = view;

%% Obtain an A given a wantedpSLocation by linear interpolation

AInterp = zeros(4,4);
A1stInterp = zeros(4,4);
A2ndInterp = zeros(4,4);
for i = 1:4
    for j = 1:4   
        coefValues = AComplete(i,j,:);
        coefValues = coefValues(:);
        yi = interp1(pSY,coefValues, wantedPSLocation);
        AInterp(i,j) = yi;
        
        coefValues = A1stComplete(i,j,:);
        coefValues = coefValues(:);
        yi = interp1(pSY,coefValues, wantedPSLocation);
        A1stInterp(i,j) = yi;
        
        coefValues = A2ndComplete(i,j,:);
        coefValues = coefValues(:);
        yi = interp1(pSY,coefValues, wantedPSLocation);
        A2ndInterp(i,j) = yi;
    end
end

%% Compute ground truth LF at the WANTED Point source Field Height, and produce PSF

%assign the full coordinates for the wanted pS location to pointSource
pointSource = pSLocations(1,:);
pointSource(2) = wantedPSLocation;
[ppsf x b bMiddle xOrig bOrig ppsfCamera] = s3dVOLTRTOnePoint(pointSource, film, lens);    
    
close all;

% Plot phase space and visual PSF of linear interpolation model output
ppsfCamera = s3dVOLTCreatePSFFromLF(ppsfCamera, b)
oi = ppsfCamera.oiCreate;
    vcAddObject(oi); oiWindow;
    plotOI(oi,'illuminance mesh log');
    
%% Calculate the same result as above (PSF), but using the INTERPOLATED A Matrix instead
% The result of this experiment should be almost the same as the one above.
%
% However, we are "cheating" in a sense because we are using the ground
% truth effective aperture information to determine which rays make it
% through the lens.

close all;
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
    plotOI(oi,'illuminance mesh log');


%% Calculate the same result as above (PSF), using the 2 interpolated A matrices instead
% Once again, this result should not be too much different from the ground
% truth PSF.  The 2 interpolated A matrix split the lens into 2 halves, the
% first one from the front-most lens element to the aperture, and the
% second one from the aperture to the back-most lens element.
%
% However, we are "cheating" in a sense because we are using the ground
% truth effective aperture information to determine which rays make it
% through the lens.

close all;

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
% ppsfCamera = s3dVOLTCreatePSFFromLF(ppsfCamera, bEstInterp)
% oi = ppsfCamera.oiCreate;
%     vcAddObject(oi); oiWindow;
%     plotOI(oi,'illuminance mesh log');

%% Calculate the PSF, using the 2 A matrices instead, and the aperture in the middle
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

%*** change this parameter to change the rendered middle aperture for the
%lens
adjustedMiddleAperture = 4;
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
    plotOI(oi,'illuminance mesh log');