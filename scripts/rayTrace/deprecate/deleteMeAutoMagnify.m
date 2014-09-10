% Ray-tracing auto magnify
% This is a 2-pass process.  The first pass gives a crude rendering of the
% PSF.  Then, the system analyzes this data and figures out a good center,
% and optionally the width and height of the 2nd pass.  The 2nd pass then
% renders the PSF under these modified conditions for a good PSF rendering
% at a good scale.  This will be very important for automating the PSF
% acquisition later. 

% AL Vistalab, 2014
%% 
s_initISET
%% point sources

%for now - assume only 1 light source for simplicity
pointSources = [ 0 0 -20000];  %large distance test

%% film properties - 
film = filmObject([24 24 36.4],[1 1], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance

%% lens properties
diffractionEnabled = false;

%initialize to default
lens = lensRealisticObject([],[],[],[], 2, 50, [], diffractionEnabled);

%read lens from file
lens.readLensFile(fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat'))

%intialize ray=tracing
lens.calculateApertureSample([51 51]);

%lens illustration
lens.drawLens();

%% loop through all point sources and render PSF
%TODO - put this into function - perhaps a "camera" object?
% vcNewGraphWin; %for illustration
for curInd = 1:size(pointSources, 1);
    %calculate the origin and direction of the rays
    %     rays.traceSourceToLens(pointSources(curInd, :), lens);
    rays = lens.rayTraceSourceToLens(pointSources(curInd, :));
    
    %duplicate the existing rays, and creates one for each
    %wavelength
    rays.expandWavelengths(film.wave);
    
    %lens intersection and raytrace
    lens.rayTraceThroughLens(rays);

    %intersect with "film" and add to film
    rays.recordOnFilm(film);
end

%% Show the image

% vcNewGraphWin;
% imshow(film.image/ max(film.image(:)));

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

vcAddAndSelectObject(oi); oiWindow;

%% calculate spread statistics to figure out the 2nd pass sensor size and center

%convert the oi image to grayscale
rgbImage = oiGet(oi, 'rgbImage');
grayImage = rgb2gray(rgbImage);
figure; imshow(grayImage);
flippedGrayImage = grayImage(size(grayImage,1):-1:1, :);
figure; imshow(flippedGrayImage);

%normalize the gray image for a weighted average
flippedGrayImage = flippedGrayImage./sum(flippedGrayImage(:));

%binary image
% binImage = uint8(grayImage > 0);
% figure; imshow(binImage);
%find the mean and standard deviation
% regionprops(binImage, 'centroid')

%calculate the weighted centroid/center-of-mass
[filmDistanceX filmDistanceY] = meshgrid(linspace(-film.size(1)/2, film.size(1)/2, film.resolution(1)),  linspace(-film.size(2)/2, film.size(2)/2, film.resolution(2)));
distanceMatrix = sqrt(filmDistanceX.^2 + filmDistanceY.^2);
centroidX = sum(sum(flippedGrayImage .* filmDistanceX))
centroidY = sum(sum(flippedGrayImage .* filmDistanceY))

%figure out a way to compute the approximate spread later
%use standard deviation?  
% stdevX = std(sum(flippedGrayImage, 2))
% stdevY = std(sum(flippedGrayImage, 1))
% 
% stdevX = std((flippedGrayImage))

%for now - the width of the new sensor is set manually - this should be
%somewhat automated in the future
newWidth = .3;

%% ------ render the PSF a 2nd time using the prescribed ----------

%% film properties 
film = filmObject([centroidX centroidY 36.4],[newWidth newWidth], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance

%% intialize ray=tracing

%use more samples this time for a high quality render
lens.calculateApertureSample([101 101]);

%% loop through all point sources and render PSF
%TODO - put this into function - perhaps a "camera" object?
% vcNewGraphWin; %for illustration
for curInd = 1:size(pointSources, 1);
    %calculate the origin and direction of the rays
    %     rays.traceSourceToLens(pointSources(curInd, :), lens);
    rays = lens.rayTraceSourceToLens(pointSources(curInd, :));
    
    %duplicate the existing rays, and creates one for each
    %wavelength
    rays.expandWavelengths(film.wave);
    
    %lens intersection and raytrace
    lens.rayTraceThroughLens(rays);

    %intersect with "film" and add to film
    rays.recordOnFilm(film);
end

%% Show the image

% vcNewGraphWin;
% imshow(film.image/ max(film.image(:)));

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

vcAddAndSelectObject(oi); oiWindow;
