%% ray-tracing for realistic lens - PPSF
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
%  This specific script renders 2 PPSFs at a set distance, but with
%  different field positions.  We attempt to interpolate a PPSF for a
%  position half-way between these 2, and check the results with the ground
%  truth.
% AL Vistalab, 2014
%%
s_initISET


%% first PPSF at the center
%% point sources
pointSourceDepth = 20000;
pointSourceFieldHeight = 0;
pointSources = [ 3000 0 -pointSourceDepth];  %large distance test
% pointSources = [ 0 0 -60];  %short distance test

%% film properties -

%36.4
film = pbrtFilmC([0 0 40 ],[40 40], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance

%% lens properties
diffractionEnabled = false;   
%turning on diffraction does NOT make sense yet since we have not modeled 
%the transformation of uncertainty from the middle aperture to the end of the lens

%initialize to default
%lensRealisticObject(elOffset, elRadius, elAperture, elN, aperture, focalLength, center, diffractionEnabled, wave)
lens = lensRealisticObject([],[],[],[], 8, 50, [], diffractionEnabled);

%read lens from file
lens.readLensFile(fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat'))

%intialize ray=tracing
randJitter = false;
lens.calculateApertureSample([301 301], randJitter);

%lens illustration
lens.drawLens();
%% ray trace and save ppsf
s_s3dRayTracePPSFInterpHelper;
% save ppsf
firstRays = ppsfC();
firstRays.makeDeepCopy(modifyRays);


%% 2nd PPSF at the center to the side
%% point sources
pointSourceDepth = 20000;
pointSourceFieldHeight = 0;
pointSources = [ 6000 0 -pointSourceDepth];  %large distance test
% pointSources = [ 0 0 -60];  %short distance test

%% film properties -
film = pbrtFilmC([0 0 40 ],[40 40], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance

%% lens properties
diffractionEnabled = false;   
%turning on diffraction does NOT make sense yet since we have not modeled 
%the transformation of uncertainty from the middle aperture to the end of the lens

%initialize to default
%lensRealisticObject(elOffset, elRadius, elAperture, elN, aperture, focalLength, center, diffractionEnabled, wave)
lens = lensRealisticObject([],[],[],[], 8, 50, [], diffractionEnabled);

%read lens from file
lens.readLensFile(fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat'))

%intialize ray=tracing
randJitter = false;
lens.calculateApertureSample([301 301], randJitter);

%lens illustration
lens.drawLens();
%%
s_s3dRayTracePPSFInterpHelper;
% save ppsf
secondRays = ppsfC();
secondRays.makeDeepCopy(modifyRays);




%% 3rd PPSF in between
%% point sources
pointSourceDepth = 20000;
pointSourceFieldHeight = 0;
pointSources = [ 4500 0 -pointSourceDepth];  %large distance test
% pointSources = [ 0 0 -60];  %short distance test

%% film properties -
film = pbrtFilmC([0 0 40 ],[40 40], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance

%% lens properties
diffractionEnabled = false;   
%turning on diffraction does NOT make sense yet since we have not modeled 
%the transformation of uncertainty from the middle aperture to the end of the lens

%initialize to default
%lensRealisticObject(elOffset, elRadius, elAperture, elN, aperture, focalLength, center, diffractionEnabled, wave)
lens = lensRealisticObject([],[],[],[], 8, 50, [], diffractionEnabled);

%read lens from file
lens.readLensFile(fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat'))

%intialize ray=tracing
randJitter = false;
lens.calculateApertureSample([301 301], randJitter);

%lens illustration
lens.drawLens();
%%
s_s3dRayTracePPSFInterpHelper;
% save ppsf
middleRays = ppsfC();
middleRays.makeDeepCopy(modifyRays);


%% Interpolate between 1st and 2nd PPSF to try to produce 3rd one
averageRays = ppsfC();
averageRays.makeDeepCopy(firstRays);
%average between the first 2 ppsfs
averageRays.origin = (firstRays.origin + secondRays.origin)./2;
averageRays.direction = (firstRays.direction + secondRays.direction)./2;
averageRays.apertureLocation = (firstRays.apertureLocation + secondRays.apertureLocation)./2;

% ray-trace the last bit - from lens to sensor
%modify the film and see the consequences on the PSF - these computations
%should be very fast
modifyRays = ppsfC();
modifyRays.makeDeepCopy(averageRays);
% 
% newRadius = 2;
outsideAperture = [];
% outsideAperture = modifyRays.apertureLocation(:,1).^2 + modifyRays.apertureLocation(:,2).^2 > newRadius^2;

%modify so only x >0 shows up
% outsideAperture = modifyRays.apertureSamples.X > 0; 

%remove outside of aperture elements
%TODO: make this into a function
modifyRays.origin(outsideAperture, : ) = [];   %this needs to be fixed later
modifyRays.direction(outsideAperture, : ) = [];
modifyRays.wavelength(outsideAperture) = [];
modifyRays.waveIndex(outsideAperture) = [];
modifyRays.apertureLocation(outsideAperture, :) = [];
modifyRays.apertureSamples.X(outsideAperture) = []; 
modifyRays.apertureSamples.Y(outsideAperture) = [];

filmCell = cell(1,1);
%first try at 36.4 sensor distance
filmCell{1} = pbrtFilmC([0 0 40 ],[40 40], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance
%intersect with "film" and add to film
disp('-----record on film-----');
tic
modifyRays.recordOnFilm(filmCell{1});
toc


% Show the images


% vcNewGraphWin;
% imshow(film.image/ max(film.image(:)));

for i = 1:length(film)

    oi = oiCreate;
    oi = initDefaultSpectrum(oi);
    oi = oiSet(oi, 'wave', filmCell{i}.wave);
    oi = oiSet(oi,'photons',filmCell{i}.image);


    optics = oiGet(oi,'optics');
    optics = opticsSet(optics,'focal length',lens.focalLength/1000);
    optics = opticsSet(optics,'fnumber', lens.focalLength/(2*1));
    oi = oiSet(oi,'optics',optics);
    hfov = rad2deg(2*atan2(filmCell{i}.size(1)/2,lens.focalLength));
    oi = oiSet(oi,'hfov', hfov);
    
    temp = filmCell{i}.position;
    filmDistance = temp(3);
    oi = oiSet(oi, 'name', ['filmDistance: ' num2str(filmDistance)]);
    vcAddAndSelectObject(oi); oiWindow;
end


