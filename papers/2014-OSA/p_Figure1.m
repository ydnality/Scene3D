%% p_Figure1 for 2014 OSA Conference
%
%  In the ray trace calculation imagine the z-axis is positive on the right
%  and negative on the left.  The scene is located on the left.  The lens
%  is placed near zero (and has some thickness).  The film is on the right
%  (z is positive).
%
%  
% AL

%%
s_initISET

%% Initialize scene, lens and film properties

% We will loop through the lens positions
% pX = 0; pY = 0; pZ = -20000;   % millimeters
pX = 0:-400:-800; pY = 0; pZ = [-10000 -8000];   % millimeters
[X, Y, Z] = meshgrid(pX,pY,pZ);
pointSources = [X(:), Y(:), Z(:)];

numDepths = length(pZ);
numFieldHeights = length(pX) * length(pY);
psfsPerDepth = length(pointSources)/numDepths;

% pointSources = [pX pY pZ];     % large distance test

% Create the film plane
wave = 400:100:700;            % Wavelength
wList = 1:length(wave);
fX = 0; fY = 0; fZ = 51.5;       % mm

% Film width and height
fW = 20;  % mm
fH = 20;  % mm

numPixelsW = 201;
numPixelsH = 201;
%for now - the width of the new sensor is set manually - this should be
%somewhat automated in the future
newWidth = .3;    %%mm

%declare film - put later for the loop

%% Describe the lens

% At some point, we will read in a lens description file

% Multicomponent lens properties
% This goes from the light through the lens to the film++++++
offset   = [0 1.5 1.5];   % Distances between surfaces (deltas) 
radius   = [67 0 -67];    % Radius of curvature, 0 means aperture
aperture = [3 2 3];       % Circular apertures, these are the radii in mm

% Index of refraction to the right of each surface
%(ray.wavelength - 550) * -.04/(300) + curEl.n;

firstN = (wave - 550) * -.04/(300) + 1.65; %linearly changes the 1.65 material
n = [firstN' zeros(4, 1) ones(4,1)]; %index of refraction (wavelength x element)

nSamples = 51;           % On the first aperture. x,y, before cropping
nSamplesHQ = 101;
% May not be needed ... AL
lX = 0; lY = 0; lZ = -1.5;
lensCenterPosition = [lX lY lZ];  % Eventually calculate this given the lens file


idx = find(radius==0);  % This is the middle of the lens aperture size
fLength = 50;           % mm.  We should derive this using the lensmaker's equation
% For multiple lenses, we add up the power using something from the web

diffractionEnabled = false;
lens = lensRealisticObject(offset,radius,aperture, n, aperture(idx), fLength, lensCenterPosition, diffractionEnabled, wave);
lens.calculateApertureSample([nSamples nSamples]);

%% Loop on wavelength and depth to create PSFs 
% These psfs will be for different field heights, depths, and wavelengths
oiList = cell(1,size(pointSources, 1))
for curInd = 1:size(pointSources, 1)
    %---initial low quality render
    film = filmObject([fX fY fZ],[fW fH], wave, [wave(:) wList(:)], [numPixelsW numPixelsH length(wave)]);
    psfCamera = psfCameraObject(lens, film, pointSources(curInd, :));
    oi = psfCamera.estimatePSF();
    vcAddObject(oi); oiWindow;
    
    %---figure out center position
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

    %---re-render
    film = filmObject([centroidX centroidY fZ],[newWidth newWidth], wave, [wave(:) wList(:)], [numPixelsW numPixelsH length(wave)]);   %large distance
    %use more samples this time for a high quality render
    lens.calculateApertureSample([nSamplesHQ nSamplesHQ]);
    psfCamera = psfCameraObject(lens, film, pointSources(curInd, :));
    oiList{curInd} = psfCamera.estimatePSF();
    vcAddObject(oiList{curInd}); oiWindow;
end

%% Show PSFs for figure

%form PSF matrix
PSF = zeros(numPixelsW, numPixelsH, length(wave), numDepths, numFieldHeights);
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

%make figure
for waveInd = wList 
    for depthInd = 1:numDepths
        PSFMosaic = [];
        for fHIndex= 1:psfsPerDepth
            curPSF = PSF(:,:, waveInd, depthInd, fHIndex);
            PSFMosaic = [PSFMosaic curPSF];  %concatenate PSFs for visualization
        end
        figure; imshow(PSFMosaic/max(PSFMosaic(:))); 
        title(['Depth:' num2str(pZ(depthInd)) '; Wavelength: ' num2str(wave(waveInd)) ]); 
    end
end
%%


%%