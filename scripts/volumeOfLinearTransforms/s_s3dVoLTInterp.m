%% s_s3dVoLTInterp
%
% Volume of Linear Transforms (VOLT) ray-tracing for multi-element lenses
%
% Examine the interpolation properties of the linear transforms
%
% AL Vistalab, 2014
%%
s_initISET

%% Specify different point source positions in the object volume

% We compute a linear transform for three field heights
pSY = [1.5 1.7 1.9];
pSZ = [-102.75];   %values must be monotonically increasing!!
wave = 550;

%% Define the Lens and Film

% lensFileName = fullfile(s3dRootPath,'data', 'lens', 'dgauss.50mm.dat');
lensFileName = fullfile(s3dRootPath,'data', 'lens', '2ElLens.dat');

nSamples = 151;
apertureMiddleD = 8;   % mm
lens = lensC('apertureSample', [nSamples nSamples], ...
    'fileName', lensFileName, ...
    'apertureMiddleD', apertureMiddleD);

% We don't want to interpolate
lens.wave = wave;

% position - relative to center of final lens surface
% size - 'mm'
% wavelength samples
film = pbrtFilmC('position', [0 0 100 ], ...
    'size', [10 10], ...
    'wave', wave);

%% Compute VOLT model
VoLTObject = VoLTC('lens', lens, 'film', film, 'fieldPositions', pSY, 'depths', pSZ, 'wave', wave); 
VoLTObject.calculateMatrices();

vcNewGraphWin([],'tall');
for ii=1:3
    subplot(3,1,ii); surf(VoLTObject.ACollection(:,:,ii));
end

0.5*(VoLTObject.ACollection(:,:,1) + VoLTObject.ACollection(:,:,3)) - ...
    VoLTObject.ACollection(:,:,2)

%%  Compare the PSF from the interpolated and ray-traced

% First, take the rays at the entrance pupil and transform them to the
% middle aperture.

lt1 = VoLTObject.A1stCollection(:,:,2);
apSize = lens.get('middle aperture d');

middleXY = ppsf.aEntranceInt.XY;
withinAperture = middleXY(:,1).^2 + middleXY(:,2).^2 <= adjustedMiddleAperture.^2;%apertureMiddleD/2;
%middleAperture = diag(middleAperture);

firstHalf = A1stInterp * xOrig;
firstHalfBlock = firstHalf(:, withinAperture); %Apply aperture to rays from the first half
bEstInterp = A2ndInterp * firstHalfBlock;

bOrigCropped = bOrig(:, withinAperture); %ground truth rays




    xCurrentWave = xOrig(:,inCurrentWaveBand);  %sifts out rays that only of the current wave band
    middleXYCurrentWave = middleXY(inCurrentWaveBand,:);
    withinAperture = middleXYCurrentWave(:,1).^2 + middleXYCurrentWave(:,2).^2 <= adjustedMiddleApertureRadius.^2;
    
    A1stInterpCurrentWave = A1stInterp(:,:,w);  %use 2 dimensions of the matrices for multiplication
    A2ndInterpCurrentWave = A2ndInterp(:,:,w);
    firstHalf = A1stInterpCurrentWave * xCurrentWave;
    firstHalfBlock = firstHalf(:, withinAperture); %Apply aperture to rays from the first half
    bEstInterp = A2ndInterpCurrentWave * firstHalfBlock;
    
    bOrigMaskedCurrentWave = bOrig(:, withinAperture); %ground truth rays
    finalB = cat(2, finalB, bEstInterp);   %concatenate interpolatedB to final B list
    waveIndex = cat(1, waveIndex, ones(size(bEstInterp,2), 1)* w); 

%% Obtain an A given a wanted pSLocation by linear interpolation
% A different A matrix will be calculated for each wavelengths.  We will
% loop through all the wavelengths.  The wave samples assumed are the ones
% from the lens.  These should be synchronized with everything else.  

%wantedWavelength = 550;

%full model from front-most element to back-most element
AInterp = zeros(4,4, length(wave));  
%half the model from the front-most element to the middle aperture
A1stInterp = zeros(4,4, length(wave));  
%the second half of the model from the middle aperture to the back-most element
A2ndInterp = zeros(4,4, length(wave));  

for w = 1:length(wave)
    [AInterp1Wave, A1stInterpCurrentWave, A2ndInterpCurrentWave ] = VoLTObject.interpolateA(wantedPSLocation, wave(w));
    AInterp(:,:,w) = AInterp1Wave;
    A1stInterp(:,:,w) = A1stInterpCurrentWave;
    A2ndInterp(:,:,w) = A2ndInterpCurrentWave;
end


%% Compute ground truth LF at the WANTED Point source Field Height

% For this moment, we only run at one depth
% We set the wanted field height to the second coordinate,
% The first coordinate is always 0.
pointSource = wantedPSLocation;

% Compute the plenoptic pointspread that has the light field information
[ppsf, x, b, bMiddle, xOrig, bOrig, ppsfCamera] = ...
    s3dVOLTRTOnePoint(pointSource, film, lens);

close all;

%% Now compute the PSF

% Plot phase space and visual PSF of linear interpolation model output
ppsfCamera = s3dVOLTCreatePSFFromLF(ppsfCamera, b);

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
adjustedMiddleApertureRadius = 4;
close all;

middleXY = ppsf.aEntranceInt.XY;
%withinAperture = middleXY(:,1).^2 + middleXY(:,2).^2 <= adjustedMiddleApertureRadius.^2;%apertureMiddleD/2;
finalB = [];
waveIndex = [];
%loops through all different waves
for w = 1:length(wave)
    
    inCurrentWaveBand = (ppsf.waveIndex == w);
    xCurrentWave = xOrig(:,inCurrentWaveBand);  %sifts out rays that only of the current wave band
    middleXYCurrentWave = middleXY(inCurrentWaveBand,:);
    withinAperture = middleXYCurrentWave(:,1).^2 + middleXYCurrentWave(:,2).^2 <= adjustedMiddleApertureRadius.^2;
    
    A1stInterpCurrentWave = A1stInterp(:,:,w);  %use 2 dimensions of the matrices for multiplication
    A2ndInterpCurrentWave = A2ndInterp(:,:,w);
    firstHalf = A1stInterpCurrentWave * xCurrentWave;
    firstHalfBlock = firstHalf(:, withinAperture); %Apply aperture to rays from the first half
    bEstInterp = A2ndInterpCurrentWave * firstHalfBlock;
    
    bOrigMaskedCurrentWave = bOrig(:, withinAperture); %ground truth rays
    finalB = cat(2, finalB, bEstInterp);   %concatenate interpolatedB to final B list
    waveIndex = cat(1, waveIndex, ones(size(bEstInterp,2), 1)* w); %keep track of waveIndex
end


% % Calculate errors
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

% Visualize PSF and phase space
ppsfCamera = s3dVOLTCreatePSFFromLF(ppsfCamera, finalB, withinAperture, waveIndex)
oi = ppsfCamera.oiCreate;
vcAddObject(oi); oiWindow;
plotOI(oi,'illuminance mesh linear');
%% End 