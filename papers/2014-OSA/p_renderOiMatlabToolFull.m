%% 2014 OSA Conference: Full Forward calculation
%
% This script takes a scene image with depth map, and applies the series of
% PSF's to it, depending on the position, wavelength, and depth.



%%
s_initISET

%% Render a pinhole scene.  It can be treated as the scene radiance

%TODO: save these scenes and depth maps so that brian can import them
sceneName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'indestructibleObject', 'mainPinhole.pbrt');
scene = s3dRenderScene(sceneName, 'indObj');
vcAddObject(scene); sceneWindow;

%% Render the depthmap for that pinhole scene
sceneName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'indestructibleObject', 'mainPinholeDepth.pbrt');
depthMap = s3dRenderDepthMap(sceneName, 1);
figure; imagesc(depthMap);
title('Depth-map For Indestructible Object(units mm)');
colorbar;
scene = sceneSet(scene, 'depthmap', depthMap);

%% If modded pbrt is NOT installed on this system, run this command to 
% load a scene file

sceneFileName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'indestructibleObject', 'pinholeSceneFile.mat');
scene = load(sceneFileName);
scene = scene.scene;
vcAddObject(scene); sceneWindow;
%% Render a collection of PSFs (similar to p_Figure1.m) - could make this a function later or separate script.
% For each PSF rendered, we will calculate the centroid of that PSF, and
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

    %% Describe the point sources 

    % We will loop through the lens positions
    pX = 0:-800:-2400; pY = 0; pZ =[-70 -80 -90 -110];% millimeters
    normalizingZ = -16000; %mm assumed reference Z point.  All other points will use the same field angle as this reference Z point
    [X, Y, Z] = meshgrid(pX,pY,pZ);
    
    %adjust for approximate difference in field position when Z changes
    for i = 1:length(pZ)
        X(:,:,i) = X(:,:,i) *  pZ(i)/normalizingZ; 
    end
    pointSources = [X(:), Y(:), Z(:)];

    numDepths = length(pZ);
    numFieldHeights = length(pX) * length(pY);
    psfsPerDepth = size(pointSources, 1)/numDepths;
    jitterFlag = true;   %enable jitter for lens front element aperture samples
    nLines = false;  %number of lines to draw for debug illustrations.  

    %%  Declare film properties (for PSF recording only, NOT for the forward calculation)
    wave = 400:100:700;            % Wavelength
    wList = 1:length(wave);
    fX = 0; fY = 0; fZ = 135.5;       % mm.  Film position
    
    % Film width and height
    fW = 60;  % mm
    fH = 60;  % mm
    % Film resolution (preview)
    numPixelsW = 151;
    numPixelsH = 151;
    % Film resolution (final render)
    numPixelsWHQ = 75;
    numPixelsHHQ = 75;

    %for now - the width of the high quality sensor is set manually - this should be
    %somewhat automated in the future
    newWidth = 10;    %%mm

    %% Describe the lens

    % Multicomponent lens properties
    % This goes from the light through the lens to the film
    zPos   = [-3 -1.5 0];   % Z intercept positions of lens surfaces
    radius   = [67 0 -67];    % Radius of curvature, 0 means aperture
    aperture = [10 10 10];       % Circular apertures, these are the radii in mm

    % Index of refraction to the right of each surface
    %(ray.wavelength - 550) * -.04/(300) + curEl.n;
    firstN = (wave - 550) * -.04/(300) + 1.65; %linearly changes the 1.65 material
    n = [firstN' zeros(length(wave), 1) ones(length(wave),1)]; %index of refraction (wavelength x element)

    nSamples = 25;           % On the first aperture. x,y, before cropping
    nSamplesHQ = 801;        % Number of samples for the HQ render
    diffractionEnabled = false;    %disable diffraction for this example
    idx = find(radius==0);  % This is the middle of the lens aperture size mm
    fLength = 50;           % Todo: We should derive this using the lensmaker's equation
    % For multiple lenses, we add up the power using something from the web

    % Populate lens surface array ussing given property arrays
    lensSurfaceArray = lensSurfaceObject();
    for i = 1:length(zPos)
        lensSurfaceArray(i) = lensSurfaceObject('sRadius', radius(i), 'apertureD', aperture(i), 'zPos', zPos(i), 'n', n(:, i));
    end

    % Declare lens
    lens = lensMEObject('surfaceArray', lensSurfaceArray, 'focalLength', fLength, 'diffractionEnabled', diffractionEnabled, 'wave', wave, 'aperturesample', [nSamples nSamples]);
    lens.apertureMiddleD = 10;
    % lens.calculateApertureSample([nSamples nSamples]);

    %% Loop on wavelength and depth to create PSFs 
    % These psfs will be for different field heights, depths, and wavelengths
    oiList = cell(1,size(pointSources, 1))
    for curInd = 1:size(pointSources, 1)
    % for curInd = 1:1
        %---initial low quality render
        film = pbrtFilmObject('position', [fX fY fZ], 'size', [fW fH], 'wave', wave, 'resolution', [numPixelsW numPixelsH length(wave)]);
        psfCamera = psfCameraObject('lens', lens, 'film', film, 'pointsource', pointSources(curInd, :));
        oi = psfCamera.estimatePSF();
        %         vcAddObject(oi); oiWindow;

        %---figure out center pos by calculating the centroid of gray image
        %convert the oi image to grayscale
        rgbImage = oiGet(oi, 'rgbImage');
        grayImage = rgb2gray(rgbImage);
        %         figure; imshow(grayImage);
        grayImage = grayImage./sum(grayImage(:));
        flippedGrayImage = grayImage(size(grayImage,1):-1:1, :);
        %         figure; imshow(flippedGrayImage);
        %normalize the gray image for a weighted average
        flippedGrayImage = flippedGrayImage./sum(flippedGrayImage(:));
        %calculate the weighted centroid/center-of-mass
        [filmDistanceX filmDistanceY] = meshgrid(linspace(-film.size(1)/2, film.size(1)/2, film.resolution(1)),  linspace(-film.size(2)/2, film.size(2)/2, film.resolution(2)));
        distanceMatrix = sqrt(filmDistanceX.^2 + filmDistanceY.^2);
        centroidX = sum(sum(flippedGrayImage .* filmDistanceX)); 
        centroidY = sum(sum(flippedGrayImage .* filmDistanceY));

        %---re-render image under new center position and width
        smallFilm = pbrtFilmObject('position', [centroidX centroidY fZ], 'size', [newWidth newWidth], 'wave', wave, 'resolution', [numPixelsWHQ numPixelsHHQ length(wave)]);
        %use more samples this time for a high quality render
        lens.apertureSample = ([nSamplesHQ nSamplesHQ]);
        psfCamera = psfCameraObject('lens', lens, 'film', smallFilm, 'pointsource', pointSources(curInd, :));
        oiList{curInd} = psfCamera.estimatePSF(nLines, jitterFlag);
%         vcAddObject(oiList{curInd}); oiWindow;
    end

    %% Compute PSF collection matrix. 
    % This will serve as a lookup table for later parts of the script

    % Form PSF matrix
    PSF = zeros(numPixelsWHQ, numPixelsHHQ, length(wave), numDepths, numFieldHeights);
    for waveInd = wList 
        for depthInd = 1:numDepths
            for fHIndex = 1:psfsPerDepth
                longInd = (depthInd - 1) * psfsPerDepth + fHIndex;
                curOi = oiList{longInd};
                curPhotons = oiGet(curOi, 'photons');
                curPSF = curPhotons(:,:, waveInd);
                PSF(:,:,waveInd,depthInd, fHIndex) = curPSF; %put PSF for current depth and wavelength in the matrix;
            end
        end
    end

    % Key data to know for interpolation later
    PSFFieldHeightSamples = atan(pX/normalizingZ) * 180/pi
    PSFDepthSamples = -Z(1,1,:);
    PSFDepthSamples = PSFDepthSamples(:)
    PSFStructure.fHAngle = PSFFieldHeightSamples;
    PSFStructure.depth = PSFDepthSamples';
    PSFStructure.wave = wave;
    PSFStructure.PSF = PSF;
    PSFStructure.film = smallFilm;



%% Apply proper PSF for every pixel - this will become a function later

% Dimensions of scene data
numRows = sceneGet(scene, 'rows');
numCols = sceneGet(scene,'cols');
% Size of PSF film pixel
smallPixelSize = PSFStructure.film.size(1)/size(PSFStructure.film.image,1);     %mm
% Size of scene sample pixel
largePixelSize = 70/sqrt(2)/numRows;  %mm   
renderWave = 400:100:700;
filmDistance = PSFStructure.film.position(3); %should match the PSFStructure film Z to be consistent


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

% Loop through wavelength
for waveInd = 1:length(renderWave)
    curWave = renderWave(waveInd);
    waveInd
    % Loop through rows and cols
    for ii = 1:numRows
        for jj = 1:numCols

            %for each pixel, 
            %-figure out the field height(in angles)
            
            %-figure out the PSF for that field position and that depth
            %-reposition the PSF at that pixel and multiply it by
            %the pixel value.
            %-add to the sum
            
            % Calculate field height (in degrees)
            hypotenuse = sqrt((ii - numRows/2)^2 + (jj - numCols)^2) * largePixelSize;
            curFieldHeightAngle = atan(hypotenuse/filmDistance) * 180/pi; 

            depth = dM(ii, jj);
            
            % Calculates the correct PSF for the current field angle, depth,
            % wavelength, given the PSf structure
            currentPSF = s3dPSFLookUp(curFieldHeightAngle, depth, curWave,PSFStructure);
            
            % Scale PSF to the right size
            scaledCPSF = imresize(currentPSF, scaleFactor);
            scaledCPSF = scaledCPSF./sum(scaledCPSF(:)); %normalize
            scaledPSFNumRows = size(scaledCPSF, 1);
            
            % Define center and starting window
            centerPos = ceil(size(scaledCPSF, 1)/2);
            
            startingRow = ii + padAmount - centerPos;
            endingRow  = startingRow + scaledPSFNumRows-1;
            startingCol = jj + padAmount - centerPos;
            endingCol = startingCol + scaledPSFNumRows-1;
            
            % Add to the sum
            photonSum(startingRow:endingRow,startingCol:endingCol, waveInd) = ...
                photonSum(startingRow:endingRow,startingCol:endingCol, waveInd)+ ...
                unBlurredPhotons(ii + padAmount,jj + padAmount, waveInd)*scaledCPSF;
        end
    end
end

% Assign sum as oi values
oi = oiSet(oi, 'wave', renderWave);
oi = oiSet(oi,'cphotons',photonSum);
vcAddObject(oi); oiWindow;

% debugging 
% figure; 
% imshow(photonSum(:,:,1)./ max(photonSum(:)) * 4)
% figure; 
% imshow(unBlurredPhotons(:,:,1)./ max(unBlurredPhotons(:)) * 4)
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
