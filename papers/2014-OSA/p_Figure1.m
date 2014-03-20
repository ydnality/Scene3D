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
pointsPerImage = length(pX);
numImages = length(pointSources)/pointsPerImage;

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

nSamples = 25;           % On the first aperture. x,y, before cropping
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

%% Create the PSF for each point

% for imageIndex = 1:numImages
    %we need a new film per image
    % TODO:  Get rid of wList part

%     curPointSources = pointSources((imageIndex-1) * pointsPerImage + 1: imageIndex * pointsPerImage, :);
    
    for curInd = 1:size(pointSources, 1);
        %---initial low quality render
        film = filmObject([fX fY fZ],[fW fH], wave, [wave(:) wList(:)], [numPixelsW numPixelsH length(wave)]);
        
        %calculate the origin and direction of the rays
        rays = lens.rayTraceSourceToLens(pointSources(curInd, :));

        %duplicate the existing rays, and creates one for each
        %wavelength
        rays.expandWavelengths(film.wave);

        %lens intersection and raytrace
        lens.rayTraceThroughLens(rays);

        %intersect with "film" and add to film
        rays.recordOnFilm(film);

        %show the oi
        oi = oiCreate;
        oi = initDefaultSpectrum(oi);
        oi = oiSet(oi, 'wave', film.wave);
        oi = oiSet(oi,'photons',film.image);
        optics = oiGet(oi,'optics');
        optics = opticsSet(optics,'focal length',lens.focalLength/1000);
        optics = opticsSet(optics,'fnumber', lens.focalLength/(2*1));
        oi = oiSet(oi,'optics',optics);
        hfov = rad2deg(2*atan2(film.size(1)/2,lens.focalLength));
        oi = oiSet(oi,'hfov', hfov);
        vcAddObject(oi); oiWindow;
        
        %---figure out center and re-render 
        
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
        centroidX = sum(sum(flippedGrayImage .* filmDistanceX)) + film.size(1)/film.resolution(1)/2; %last one is correction for rounding - see if we can remove this
        centroidY = sum(sum(flippedGrayImage .* filmDistanceY)) - film.size(1)/film.resolution(1)/2;
        
        
        %---re-render
        film = filmObject([centroidX centroidY fZ],[newWidth newWidth], wave, [wave(:) wList(:)], [numPixelsW numPixelsH length(wave)]);   %large distance
        %use more samples this time for a high quality render
        lens.calculateApertureSample([nSamplesHQ nSamplesHQ]);
        
        %calculate the origin and direction of the rays
        rays = lens.rayTraceSourceToLens(pointSources(curInd, :));
        
        %duplicate the existing rays, and creates one for each
        %wavelength
        rays.expandWavelengths(film.wave);
        
        %lens intersection and raytrace
        lens.rayTraceThroughLens(rays);
        
        %intersect with "film" and add to film
        rays.recordOnFilm(film);
        
        %show the oi
        oi = oiCreate;
        oi = initDefaultSpectrum(oi);
        oi = oiSet(oi, 'wave', film.wave);
        oi = oiSet(oi,'photons',film.image);
        optics = oiGet(oi,'optics');
        optics = opticsSet(optics,'focal length',lens.focalLength/1000);
        optics = opticsSet(optics,'fnumber', lens.focalLength/(2*1));
        oi = oiSet(oi,'optics',optics);
        hfov = rad2deg(2*atan2(film.size(1)/2,lens.focalLength));
        oi = oiSet(oi,'hfov', hfov);
        vcAddObject(oi); oiWindow;
    end
    
    

% end


%% Loop on wavelength and depth to create PSFs 

% These psfs will be for different field heights, depths, and wavelengths


%% Show pictures


%%


%%