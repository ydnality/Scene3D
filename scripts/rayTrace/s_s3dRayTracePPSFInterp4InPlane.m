%% Experimenting with plenoptic PSF (pPSF)
%
%  This uses Snell's law and a lens prescription to create a ray trace of
%  four principal points and then to try to use the pPSFs of these points
%  to interpolate the pPSF of a fifth point within the region spanned by
%  the first four.
%
%  We split the calculation of the PSF into 2 steps.
%
%    1. The first step traces the rays from a single point in the scene
%    through the lens, to the basck surface of the multiple lens optics.
%    At this point, the rays may be saved as data.
%
%    2. Next, the rays are traced from the back surface of the final lens
%    to the sensor (this process is reasonably efficient and doesn't take
%    much time). Using the pPSF format, the PSF at the sensor surface for
%    different sensor depths may be calculated.
%
%  PROGRAMMING:
%  The script is too long, and we need to start writing functions so that
%  the length is shortened and the clarity increased.
%  Wavelength is used, but no chromatic aberration yet.
%
%
%
% AL Vistalab, 2014
%%
s_initISET


%% --Specify PPSF Point Locations---
    pointSourceDepth = 1000;
    lowerRightPosition = [40 40 -pointSourceDepth];   %lower right
    lowerLeftPosition =  [50 40 -pointSourceDepth];   %lower left
    upperRightPosition = [40 50 -pointSourceDepth];   %upper right
    upperLeftPosition = [50 50 -pointSourceDepth];   %upper left
    
    distance = 10;
    newPosition = [42.5 42.5 -pointSourceDepth];
    lowerRightWeight = (distance - (newPosition(1) - lowerRightPosition(1)))/distance * (distance - (newPosition(2) - lowerRightPosition(2)))/distance
    lowerLeftWeight = (distance - -(newPosition(1) - lowerLeftPosition(1)))/distance * (distance - (newPosition(2) - lowerLeftPosition(2)))/distance
    upperRightWeight = (distance - (newPosition(1) - upperRightPosition(1)))/distance * (distance - -(newPosition(2) - upperRightPosition(2)))/distance
    upperLeftWeight = (distance - -(newPosition(1) - upperLeftPosition(1)))/distance * (distance - -(newPosition(2) - upperLeftPosition(2)))/distance
    
%% --first PPSF lower right--
    %% point sources
    pointSources = lowerRightPosition;  %large distance test
    % pointSources = [ 0 0 -60];  %short distance test

    %% film properties
    %36.4
    film = pbrtFilmObject([0 0 40 ],[10 10], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance

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
    ppsfCamera = ppsfCameraObject(lens, film, pointSources);
    ppsfRays = ppsfCamera.estimatePPSF();

    %modify the rays for any aperture changes here
    modifyRays = ppsfObject();
    modifyRays.makeDeepCopy(ppsfRays);

    %trace from end of lens to sensor
    modifyRays.recordOnFilm(ppsfCamera.film);

    %show image
    ppsfCamera.showFilm();

    % save ppsf
    lowerRightRays = ppsfObject();
    lowerRightRays.makeDeepCopy(modifyRays);

%% --2nd PPSF lower left --
    %% point sources
    pointSources = lowerLeftPosition;  %large distance test
    % pointSources = [ 0 0 -60];  %short distance test

    %% film properties -
    film = pbrtFilmObject([0 0 40 ],[10 10], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance

    %% lens properties
    diffractionEnabled = false;   
    %turning on diffraction does NOT make sense yet since we have not modeled 
    %the transformation of uncertainty from the middle aperture to the end of the lens

    %initialize to default
    lens = lensRealisticObject([],[],[],[], 8, 50, [], diffractionEnabled);

    %read lens from file
    lens.readLensFile(fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat'))

    %intialize ray=tracing
    randJitter = false;
    lens.calculateApertureSample([301 301], randJitter);

    %lens illustration
    lens.drawLens();
    %% ray trace and save ppsf
    ppsfCamera = ppsfCameraObject(lens, film, pointSources);
    ppsfRays = ppsfCamera.estimatePPSF();

    %modify the rays for any aperture changes here
    modifyRays = ppsfObject();
    modifyRays.makeDeepCopy(ppsfRays);

    %trace from end of lens to sensor
    modifyRays.recordOnFilm(ppsfCamera.film);

    %show image
    ppsfCamera.showFilm();

    % save ppsf
    lowerLeftRays = ppsfObject();
    lowerLeftRays.makeDeepCopy(modifyRays);

%% --3rd PPSF upper right--
    %% point sources
    pointSources = upperRightPosition;  %large distance test
    % pointSources = [ 0 0 -60];  %short distance test

    %% film properties -
    film = pbrtFilmObject([0 0 40 ],[10 10], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance

    %% lens properties
    diffractionEnabled = false;   
    %turning on diffraction does NOT make sense yet since we have not modeled 
    %the transformation of uncertainty from the middle aperture to the end of the lens

    %initialize to default
    lens = lensRealisticObject([],[],[],[], 8, 50, [], diffractionEnabled);

    %read lens from file
    lens.readLensFile(fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat'))

    %intialize ray=tracing
    randJitter = false;
    lens.calculateApertureSample([301 301], randJitter);

    %lens illustration
    lens.drawLens();
    %% ray trace and save ppsf
    ppsfCamera = ppsfCameraObject(lens, film, pointSources);
    ppsfRays = ppsfCamera.estimatePPSF();

    %modify the rays for any aperture changes here
    modifyRays = ppsfObject();
    modifyRays.makeDeepCopy(ppsfRays);

    %trace from end of lens tmato sensor
    modifyRays.recordOnFilm(ppsfCamera.film);

    %show image
    ppsfCamera.showFilm();

    % save ppsf
    upperRightRays = ppsfObject();
    upperRightRays.makeDeepCopy(modifyRays);
    
%% --4th PPSF upper left--
    %% point sources
    pointSources = upperLeftPosition;  %large distance test
    % pointSources = [ 0 0 -60];  %short distance test

    %% film properties -
    film = pbrtFilmObject([0 0 40 ],[10 10], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance

    %% lens properties
    diffractionEnabled = false;   
    %turning on diffraction does NOT make sense yet since we have not modeled 
    %the transformation of uncertainty from the middle aperture to the end of the lens

    %initialize to default
    lens = lensRealisticObject([],[],[],[], 8, 50, [], diffractionEnabled);

    %read lens from file
    lens.readLensFile(fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat'))

    %intialize ray=tracing
    randJitter = false;
    lens.calculateApertureSample([301 301], randJitter);

    %lens illustration
    lens.drawLens();
    %% ray trace and save ppsf
    ppsfCamera = ppsfCameraObject(lens, film, pointSources);
    ppsfRays = ppsfCamera.estimatePPSF();

    %modify the rays for any aperture changes here
    modifyRays = ppsfObject();
    modifyRays.makeDeepCopy(ppsfRays);

    %trace from end of lens tmato sensor
    modifyRays.recordOnFilm(ppsfCamera.film);

    %show image
    ppsfCamera.showFilm();

    % save ppsf
    upperLeftRays = ppsfObject();
    upperLeftRays.makeDeepCopy(modifyRays);
%% --5th PPSF middle--
    %% point sources
    pointSources = newPosition;  %large distance test
    % pointSources = [ 0 0 -60];  %short distance test

    %% film properties -
    film = pbrtFilmObject([0 0 40 ],[10 10], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance

    %% lens properties
    diffractionEnabled = false;   
    %turning on diffraction does NOT make sense yet since we have not modeled 
    %the transformation of uncertainty from the middle aperture to the end of the lens

    %initialize to default
    lens = lensRealisticObject([],[],[],[], 8, 50, [], diffractionEnabled);

    %read lens from file
    lens.readLensFile(fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat'))

    %intialize ray=tracing
    randJitter = false;
    lens.calculateApertureSample([301 301], randJitter);

    %lens illustration
    lens.drawLens();
    %% ray trace and save ppsf
    ppsfCamera = ppsfCameraObject(lens, film, pointSources);
    ppsfRays = ppsfCamera.estimatePPSF();

    %modify the rays for any aperture changes here
    modifyRays = ppsfObject();
    modifyRays.makeDeepCopy(ppsfRays);

    %trace from end of lens tmato sensor
    modifyRays.recordOnFilm(ppsfCamera.film);

    %show image
    ppsfCamera.showFilm();

    % save ppsf
    middleRays = ppsfObject();
    middleRays.makeDeepCopy(modifyRays);
    
%% --Interpolate between 4 PPSFs to try to produce 5th one--
    %% Average between 4 sets of rays
    averageRays = ppsfObject();
    averageRays.makeDeepCopy(lowerRightRays);
    %average between the first 2 ppsfs
    averageRays.origin = (lowerRightWeight * lowerRightRays.origin + lowerLeftWeight* lowerLeftRays.origin + upperLeftWeight * upperLeftRays.origin + upperRightWeight * upperRightRays.origin);
    averageRays.direction = (lowerRightWeight * lowerRightRays.direction + lowerLeftWeight* lowerLeftRays.direction + upperLeftWeight* upperLeftRays.direction + upperRightWeight*upperRightRays.direction);
    averageRays.apertureLocation = (lowerRightWeight * lowerRightRays.apertureLocation + lowerLeftWeight * lowerLeftRays.apertureLocation + upperLeftWeight * upperLeftRays.apertureLocation + upperRightWeight* upperRightRays.apertureLocation);

    % ray-trace the last bit - from lens to sensor
    %% modify the film and see the consequences on the PSF - these computations
    %should be very fast
    modifyRays = ppsfObject();
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
    
    %% record on film
    filmCell = cell(1,1);
    %first try at 36.4 sensor distance
    filmCell{1} = pbrtFilmObject([0 0 40 ],[10 10], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance
    %intersect with "film" and add to film
    disp('-----record on film-----');
    tic
    modifyRays.recordOnFilm(filmCell{1});
    toc

    %% Show the images
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


