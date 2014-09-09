%% ray-tracing for realistic lens - PPSF
%
% Towards a volume of linear transforms (volt)
%
% Check why we have a discontinuity at (0,0).  Something about the 0 always
% mapping to 0, or something ...
%
% Lens element positions are all negative, and the final exit plan can be
% considered as z = 0
%
%  This uses Snell's law and a lens prescription to create a tray trace.
%  The script is too long, and we need to start writing functions so that
%  the length is shortened and the clarity increased.
%  We are only ray-tracing ideal point sources in order to extract out point
%  spread functions.
%
%  We are also experimenting with the plenoptic point spread function
%  (PPSF).  this is a very early experiment, where we are splitting the
%  calculation into 2 steps.  The first step traces the rays from a single
%  point in the scene towards the lens.  Ray-tracing is performed through
%  the lens, and out the back aperture.  At this point, the rays may be
%  saved as data.  Next, the rays are traced from the end of the lens to
%  the sensor (this process is reasonably efficient and doesn't take much
%  time).  Using this format, differrent sensor depths may be used to
%  access the PSF.
%
%  This specific script renders 2 PPSFs at a set field position , but with
%  different depths.  We attempt to interpolate a PPSF for a
%  position half-way between these 2, and check the results with the ground
%  truth.
%
% AL Vistalab, 2014
%%
s_initISET

%% specify different point source positions
pSY = 0.01:.3:2;
pSZ = -102 * ones(length(pSY), 1);
pSLocations = [zeros(length(pSY), 1) pSY' pSZ];

%desired pSLocation for interpolation
wantedPSLocation = 1.7;


%% lens properties
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

% lens.draw();

%% compute VOLT model

[AComplete A1stComplete A2ndComplete] = s3dVOLTCreateModel(lens, film, pSLocations);

%% Interpret A matrices

% Let's interpet A. Does it agree with MP's predict calculation from the
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

%% Obtain an A given a wantedpSLocation by interpolation

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
    
%% Calculate the same result as above, but using the INTERPOLATED A Matrix instead
close all;
bEstInterp = AInterp * x;

% calculate errors
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


%% Calculate the same result as above, using the 2 A matrices instead
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

%% Calculate the same result, using the 2 A matrices instead, and the aperture in the middle

%*** change this parameter to change the rendered middle aperture for the
%lens
adjustedMiddleAperture = 4;
close all;

middleXY = ppsf.aEntranceInt.XY;
withinAperture = middleXY(:,1).^2 + middleXY(:,2).^2 <= adjustedMiddleAperture.^2;%apertureMiddleD/2;
%middleAperture = diag(middleAperture);

firstHalf = A1stInterp * xOrig;
firstHalfBlock = firstHalf(:, withinAperture);
bEstInterp = A2ndInterp * firstHalfBlock;

bOrigCropped = bOrig(:, withinAperture);

% calculate errors
% Scatter plot of positions
for ii=1:4
    ii
    vcNewGraphWin; plot(bOrigCropped(ii,:),bEstInterp(ii,:),'o');
    grid on;

    meanAbsError = mean(abs(bEstInterp(ii,:) - bOrigCropped(ii,:)));
    averageAmp = mean(abs(bOrigCropped(ii,:)));
    meanPercentErrorSplit = meanAbsError/averageAmp * 100
end

% visualize PSF and phase space
ppsfCamera = s3dVOLTCreatePSFFromLF(ppsfCamera, bEstInterp, withinAperture)
oi = ppsfCamera.oiCreate;
    vcAddObject(oi); oiWindow;
    plotOI(oi,'illuminance mesh log');

%% Future development for modifying the rays.


% Make a second ppsf object
%     modifyRays = ppsfC();
%
%     % Take the ppsfRays from the first object, copy the properties of the
%     % ppsfRays into real data, not just a pointer to the data.
%     modifyRays.makeDeepCopy(ppsfCamera.ppsfRays);
%
%     % Trace the rays lens to sensor
%     modifyRays.recordOnFilm(ppsfCamera.film);

