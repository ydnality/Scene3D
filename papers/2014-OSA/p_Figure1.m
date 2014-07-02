%% 2014 OSA Conference
%
%  This script renders the PSFs in the block diagram in Figure 1 using a 3
%  element lens (including a middle aperture).
% 
%  In the ray trace calculation imagine the z-axis is positive on the right
%  and negative on the left.  The scene is located on the left.  The lens
%  is placed near zero (and has some thickness).  The film is on the right
%  (z is positive).
%
%  
% AL

%% Initialize ISET
s_initISET

%% Describe the point sources 

% We will loop through the lens positions
% pX = 0; pY = 0; pZ = -20000;   % millimeters
pX = 0:-400:-800; pY = 0; pZ =[-16000 -8000 -4000 -2000 ];% millimeters
% pX = 0; pY = 0; pZ = [-8000];   % millimeters
[X, Y, Z] = meshgrid(pX,pY,pZ);
%adjust for approximate difference in field position when Z changes
for i = 2:length(pZ)
    X(:,:,i) = X(:,:,i) *  pZ(i)/pZ(1); 
end
pointSources = [X(:), Y(:), Z(:)];

numDepths = length(pZ);
numFieldHeights = length(pX) * length(pY);
psfsPerDepth = size(pointSources, 1)/numDepths;
jitterFlag = true;   %enable jitter for lens front element aperture samples
nLines = false;  %number of lines to draw for debug illustrations.  

%%  Declare film properties
wave = 400:100:700;            % Wavelength
wList = 1:length(wave);
fX = 0; fY = 0; fZ = 51.8145;       % mm
% Film width and height
fW = 20;  % mm
fH = 20;  % mm
% Film resolution (preview)
numPixelsW = 151;
numPixelsH = 151;
% Film resolution (final render)
numPixelsWHQ = 75;
numPixelsHHQ = 75;

%for now - the width of the new sensor is set manually - this should be
%somewhat automated in the future
newWidth = .3;    %%mm

%% Describe the lens

% Multicomponent lens properties
% This goes from the light through the lens to the film
zPos   = [-3 -1.5 0];   % Z intercept positions of lens surfaces
radius   = [67 0 -67];    % Radius of curvature, 0 means aperture
aperture = [6 4 6];       % Circular apertures, these are the radii in mm

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
lens.apertureMiddleD = 4;
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
    vcAddObject(oi); oiWindow;
    
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
    film = pbrtFilmObject('position', [centroidX centroidY fZ], 'size', [newWidth newWidth], 'wave', wave, 'resolution', [numPixelsWHQ numPixelsHHQ length(wave)]);
    %use more samples this time for a high quality render
    lens.apertureSample = ([nSamplesHQ nSamplesHQ]);
    psfCamera = psfCameraObject('lens', lens, 'film', film, 'pointsource', pointSources(curInd, :));
    oiList{curInd} = psfCamera.estimatePSF(nLines, jitterFlag);
    vcAddObject(oiList{curInd}); oiWindow;
end

%% Show PSFs for figure

%form PSF matrix
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

%key data to know for interpolation later
PSFFieldHeightSamples = atan(pX/pZ(1)) * 180/pi
PSFDepthSamples = -Z(1,1,:);
PSFDepthSamples = PSFDepthSamples(:)
PSFStructure.fHAngle = PSFFieldHeightSamples;
PSFStructure.depth = PSFDepthSamples;


% make figure
%Different field heights on 1 plot. Each depth will have an image stack of 
%different wavelengths.
colorCode = [1 0 1;   %this determines the color of the slice images
             0 0 1;
             0 1 0
             1 0 0];
         
for waveInd = wList 
    PSFMosaic = [];
    for depthInd = 1:numDepths
        rowMosaic = [];
        for fHIndex= 1:psfsPerDepth
            curPSF = PSF(:,:, waveInd, depthInd, fHIndex);
            rowMosaic = [rowMosaic curPSF];  %concatenate PSFs for visualization
        end
        PSFMosaic = rowMosaic/(max(rowMosaic(:)) * .75);
        PSFMosaic = repmat(PSFMosaic, [1 1 3]) .* repmat(reshape(colorCode(waveInd,:), [1 1 3]), size(PSFMosaic));
        
        figure; imshow(PSFMosaic);
        description = ['Depth:' num2str(pZ(depthInd)) ';Wavelength:' num2str(wave(waveInd)) ];
%         imwrite(PSFMosaic, fullfile(s3dRootPath, 'papers', '2014-OSA', [description '.png']));    
        title(description); 
    end
end

%% Old figure where field height and depth were on 1 image
for waveInd = wList 
    PSFMosaic = [];
    for depthInd = 1:numDepths
        rowMosaic = [];
        for fHIndex= 1:psfsPerDepth
            curPSF = PSF(:,:, waveInd, depthInd, fHIndex);
            rowMosaic = [rowMosaic curPSF];  %concatenate PSFs for visualization
        end
        PSFMosaic = [PSFMosaic; rowMosaic];
    end
    
    PSFMosaic = PSFMosaic/max(PSFMosaic(:));
    figure; imshow(PSFMosaic);
    
%     test = hdrRender(PSFMosaic);
%     figure; imshow(test)
    description = ['Depth:' num2str(pZ) '_Wavelength:' num2str(wave(waveInd)) ];
    imwrite(PSFMosaic, fullfile(s3dRootPath, 'papers', '2014-OSA', [description '.png']));    
    title(description); 
end