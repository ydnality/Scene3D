%% 2014 OSA Conference: Full Forward calculation
%
% This script takes a scene image with depth map, and applies the series of
% PSF's to it, depending on the position, wavelength, and depth.
%
%
% AL Copyright Vistasoft Team 2014

%%
s_initISET

%% Basic point source properties (distances)

wave = 400:100:700;            % Wavelength

% We will loop through the point positions
pX = 0:-1500:-4500;
pY = 0;
pZ =[-70 -80 -90 -100 -110];  % millimeters
normalizingZ = -16000;        % mm assumed reference Z point.
% All other points will use the same field angle as this reference Z point
        

%% Load a scene file
%
% These files were created using PBRT.  See the script s_sc3dSceneRadiance
% for an explanation of how these are created.

% Two possible scenes

% Much nicer
% sceneFileName = fullfile(s3dRootPath, 'data', 'isetScenes', 'metronome.mat');
% load(sceneFileName,'scene');

%  Much smaller
% sceneFileName = fullfile(s3dRootPath, 'data', 'isetScenes', 'textureSquare.mat');
% load(sceneFileName,'scene');

% Synthetic test and very small
% scene = sceneCreate('slanted bar',[64 64]);
% d = sceneGet(scene,'depth map');
% sz = size(d);
% d(:) = -pZ(end);
% d(1:(sz(1)/2),:) = -pZ(1);
% scene = sceneSet(scene,'depth map',d);

scene = sceneCreate('sweep frequency',96,7);
d = sceneGet(scene,'depth map');
sz = size(d);
d(:) = -pZ(end);
d(1:(sz(1)/2),:) = abs(pZ(1));
scene = sceneSet(scene,'depth map',d);

% Check scene
vcAddObject(scene); sceneWindow;

%% Render a collection of PSFs 
% Also used in p_Figure1.m - could make this a function later or separate
% script. 

% For each PSF rendered, we calculate the centroid of that PSF, and
% zoom in and calculate a "high quality" PSF using a much smaller film
% centered around that centroid.  This procedure allows us to calculate
% much higher quality PSFs, in case we wish to render very high resolution
% images.
% 
% These PSFs will be added into a PSF Collection matrix, which contains the
% whole series of PSFs depending on field height, wavelength, depth.
% 
% This collection will then be used to produce the forward calculation
% rendered image, which will use this collection as a look-up  table for
% the proper PSF to use.
%
% Please keep these indents, it helps divide the code
% hierarchically.

%%
    nWave = length(wave);
    nFH = length(pX) * length(pY);
    nDepth = length(pZ);
    nPoints = nFH*nDepth;
    
    % Old way of making point sources
    % Also, adjust for approximate difference in field position when Z changes
    %     [X, Y, Z] = meshgrid(pX,pY,pZ);
    %     for i = 1:length(pZ)
    %         X(:,:,i) = X(:,:,i) *  pZ(i)/normalizingZ;
    %     end
    
    pointSources = cell(nFH,nDepth);
    for ii=1:nFH
        for dd = 1:nDepth
            pointSources{ii,dd} = [pX(ii)*  pZ(ii)/normalizingZ, pY, pZ(dd)];
        end
    end
    
    jitterFlag   = true;   %enable jitter for lens front element aperture samples
    nLines       = false;  %number of lines to draw for debug illustrations.  

    % Plot the points in 3Space
    
    %%  Declare film properties for PSF recording (NOT for the forward calculation)

    fX = 0; fY = 0; fZ = 135.5;    % mm.  Film position
    
    % Film width and height
    fW = 120;  % mm
    fH = 120;  % mm
    
    % Film resolution (preview, large film size)
    lowRes = 151;
    
    % Film resolution (final render, small film size)
    highRes = 75;

    %for now - the width of the high quality sensor is set manually - this should be
    %somewhat automated in the future
    newWidth = 10;    %mm
    
    %% Create the lens from a lens file
    % see TwoElLens for how this lens was created
    
    lensFile = fullfile(s3dRootPath, 'data', 'lens', '2ElLens');
    load(lensFile,'lens')
    lens.apertureMiddleD = 5;

    % Preview and high quality sampling on first lens aperture
    nSamples = 25;           % On the first aperture. x,y, before cropping
    nSamplesHQ = 801;        % Number of samples for the HQ render
    
    % lens.draw;
    
    %% Loop on wavelength and depth to create PSFs
    
    % These psfs will be for different field heights, depths, and
    % wavelengths 
    wbar = waitbar(0,sprintf('Creating %i point spreads ...',numel(pointSources)));
    oiList = cell(nFH,nDepth);
    for ii = 1:nFH
        waitbar(ii/nFH,wbar);
        
        for dd = 1:nDepth

        %--- Initial low quality render
        film = pbrtFilmObject('position', [fX fY fZ], ...
            'size', [fW fH], ...
            'wave', wave, ...
            'resolution', [lowRes lowRes length(wave)]);
        
        lens.apertureSample = ([nSamples nSamples]);
        psfCamera = psfCameraObject('lens', lens, ...
            'film', film, ...
            'pointsource', pointSources{ii,dd});
        
        % What happens to each of the wavelengths?
        oi = psfCamera.estimatePSF();
        % To calculate and show, use this:
        %   oi = psfCamera.showFilm;
        %   oiGet(oi,'spatial resolution','mm')
        
        % Figure out center pos by calculating the centroid of illuminance image
        img = oiGet(oi,'illuminance');
        
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
        
        sz = oiGet(oi,'size'); mid = round(sz(1)/2);
        
        % Render image using new center position and width and higher resolution
        film = pbrtFilmObject('position', [centroidX centroidY fZ], ...
            'size', [newWidth newWidth], ...
            'wave', wave, ...
            'resolution', [highRes highRes length(wave)]);
        
        % Use more samples in the lens aperture to produce a high quality psf.
        % NOTE:  Changing the number of samples also changes the oi size.
        % This isn't good.  We need to change the sampling density without
        % changing the size.
        lens.apertureSample = ([nSamplesHQ nSamplesHQ]);
        psfCamera = psfCameraObject('lens', lens, ...
            'film', film, ...
            'pointsource', pointSources{ii,dd});
        oiList{ii,dd} = psfCamera.estimatePSF(nLines, jitterFlag);
        
        % vcAddObject(oiList{1,1}); oiWindow;
        end
    end
    
    delete(wbar)

    %% Compute PSF collection matrix. 
    % This will serve as a lookup table for later parts of the script

    % Form PSF matrix
    PSF = zeros(highRes, highRes, nWave, nDepth, nFH);
    for ww = 1:nWave 
        for dd = 1:nDepth
            for ii = 1:nFH
                PSF(:,:,ww,dd,ii) = oiGet(oiList{ii,dd}, 'photons',wave(ww));
            end
        end
    end

    % Key data to know for interpolation later.
    % This should become the ray trace structure in ISET.
    
    % PSFFieldHeightSamples = atan(pX/normalizingZ) * 180/pi;
    % PSFDepthSamples = -pZ(:);
    PSFStructure.fHAngle = atan(pX/normalizingZ) * 180/pi;
    PSFStructure.depth   = -pZ(:);
    PSFStructure.wave    = wave;
    PSFStructure.PSF     = PSF;
    PSFStructure.film    = film;
    
    % Make some plots of this PSF structure data.
    
    
    % What we really want is just this:  To make the PSFStructure here
    % compatibility with the ray trace point spread functions in ISET, and
    % then to be able to run oiCompute in raytrace mode with the depth
    % image.
    %
    % Doing this involves eliminating the code Andy copied from ISET, but
    % upgrading the ISET code to 3D.

%% This section appears to determine the size of the output image

% Dimensions of scene data
numRows = sceneGet(scene, 'rows');
numCols = sceneGet(scene,'cols');

% Size of PSF film pixel
smallPixelSize = PSFStructure.film.size(1)/size(PSFStructure.film.image,1);     %mm

% Size of scene sample pixel
largePixelSize = 70/sqrt(2)/numRows;  %mm   
renderWave = 400:100:700;

%should match the PSFStructure film Z to be consistent
filmDistance = PSFStructure.film.position(3); 

%70/sqrt(2) represents row size/diagonal  TODO: find a way to automate this
%This size should correspond to the size sensor you wish to use for the
%forward calculation.  

% Amount that PSF film needs to be scaled to be equivalent to scene sample
%size
scaleFactor = smallPixelSize/largePixelSize;

% Create an oi
oi = oiCreate;
oi = initDefaultSpectrum(oi);

% Use scene data as photons first
scene = sceneSet(scene,'wave', renderWave);
dM = sceneGet(scene, 'depthmap');
photons = sceneGet(scene, 'photons');
oi = oiSet(oi, 'wave', renderWave);
oi = oiSet(oi, 'cphotons', photons);

% Pad the oi for future processing
padAmount = (round(size(PSFStructure.film.image,1) * scaleFactor) * 2);
oi = oiPad(oi, [padAmount padAmount]);

% Obtain unblurred image, to be used later
unBlurredPhotons = oiGet(oi, 'photons');
photonSum = zeros(size(oiGet(oi,'photons')));

vcAddObject(oi); oiWindow;

%% Apply proper PSF for every pixel - this will become a function

% For each pixel,
% - Find field height, angle, wavelength, and depth
% - Find the corresponding PSF 
% - Multiply the PSF by the pixel radiance 
% - Add the result to the full image 

[X,Y] = meshgrid(1:numCols,1:numRows);
X = X - (numCols/2);
Y = Y - (numRows/2);
hyp   = sqrt( Y.^2 + X.^2)*largePixelSize;
% vcNewGraphWin; mesh(hyp);

% I think the X-axis is supposed to be zero deg.
% This isn't quite there yet.  But it is close.
ang = atan2d(Y, X);
% vcNewGraphWin; mesh(ang);

% Loop through wavelength
wBar = waitbar(0,'Slow forward rendering');
for waveInd = 1:length(renderWave)
    
    curWave = renderWave(waveInd);
    fprintf('Wavelength %.0f\n',curWave);
    
    % Loop through pixel rows and cols
    for ii = 1:numRows
        waitbar(ii/numRows,wBar);
        for jj = 1:numCols
           
            % Calculate field height (in degrees) and angle(relative to
            % center)            
            % hypotenuse = sqrt((ii - numRows/2)^2 + (jj - numCols/2)^2) * largePixelSize;
            hypotenuse = hyp(ii,jj);
            curFieldHeightAngle = atan(hypotenuse/filmDistance) * 180/pi; 
            
            % Pixel angle
            angle = ang(ii,jj);
            
            % Pixel depth
            depth = dM(ii, jj);
            
            % Calculate the PSF for the current field height, angle, depth,
            % wavelength, given the PSf structure
            currentPSF = s3dPSFLookUp(curFieldHeightAngle, depth, curWave, PSFStructure);
            
            %rotatePSF to correct orientation
            scaledCPSF = imrotate(currentPSF, angle, 'bilinear', 'crop' );
            
            % Scale PSF to the right size
            scaledCPSF = imresize(scaledCPSF, scaleFactor);
            scaledCPSF = scaledCPSF./sum(scaledCPSF(:)); %normalize

            scaledPSFNumRows = size(scaledCPSF, 1);
            scaledPSFNumCols = size(scaledCPSF, 2);
            
            % Define center and starting window
            centerPos = ceil(size(scaledCPSF, 1)/2);
            
            startingRow = ii + padAmount - centerPos;
            endingRow   = startingRow + scaledPSFNumRows - 1;
            
            startingCol = jj + padAmount - centerPos;
            endingCol   = startingCol + scaledPSFNumCols - 1;
            
            % Add to the final image 
            photonSum(startingRow:endingRow,startingCol:endingCol, waveInd) = ...
                photonSum(startingRow:endingRow,startingCol:endingCol, waveInd) + ...
                unBlurredPhotons(ii + padAmount,jj + padAmount, waveInd)*scaledCPSF;
        end
    end
end
delete(wBar);

% Assign sum as oi values
oi = oiSet(oi, 'wave', renderWave);
oi = oiSet(oi,'cphotons',photonSum);
vcAddObject(oi); oiWindow;

% % debugging 
% vcNewGraphWin; 
% imshow(photonSum(:,:,1)./ max(photonSum(:)))
% vcNewGraphWin; 
% imshow(unBlurredPhotons(:,:,1)./ max(unBlurredPhotons(:)))

%% Debugging s3dLookUp: Key data: PSF and PSFFieldHeightSamples

% wavelength = 700;
% fieldHeightAngle = 0;
% depth = 4000;
% PSFStructure
% 
% %TODO: keep track of differences in size maybe?  For now, assume app PSFs
% %are plotted on the same size sensor
% 
% currentPSF = s3dPSFLookUp(fieldHeightAngle, depth, wavelength,PSFStructure);
% vcNewGraphWin; imagesc(currentPSF); title(['PSF Depth: ' int2str(depth)]);


%% For comparison this is the backwards calculation

% If PBRT is installed, we can run this to cmpare with the backward
% calculation.

% sceneName = fullfile(s3dRootPath, 'data', 'pbrtScenes', 'metronome', 'mainDefocused2El.pbrt');
% oi = s3dRenderOI(sceneName, .050, 'indObj');
% 
% %reduce to 400:100:700
% scene = sceneCreate;
% scene = sceneSet(scene, 'photons', oiGet(oi,' photons'));
% scene = sceneSet(scene, 'wave', 400:100:700);
% 
% oi = oiSet(oi, 'wave', 400:100:700);
% oi = oiSet(oi, 'photons', sceneGet(scene, 'photons'));
% vcAddObject(oi); oiWindow;

%%
