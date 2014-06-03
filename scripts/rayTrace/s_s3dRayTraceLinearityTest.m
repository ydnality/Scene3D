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
%  This specific script renders 2 PPSFs at a set field position , but with
%  different depths.  We attempt to interpolate a PPSF for a
%  position half-way between these 2, and check the results with the ground
%  truth.
% AL Vistalab, 2014
%%
s_initISET

%% --first PPSF at the center--
    %% point sources
    pointSourceDepth = 1000;
    pointSourceFieldHeight = 0;
    pointSources = [ 0 0 -pointSourceDepth];  %large distance test
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
    ne
    ppsfRays = ppsfCamera.estimatePPSF();
    
    %project onto a plane
    ppsfRays.projectOnPlane(0);
    
    %calculate initial light field
    initialOrigin = ppsfRays.apertureSamples;
    initialDirection = zeros(size(initialOrigin.X), 3);
    initialDirection(:,1) = initialOrigin.X;
    initialDirection(:,2) = initialOrigin.Y;
    initialDirection(:,3) = -pointSourceDepth;
    initialDirection = initialDirection./repmat(sum(initialDirection,2), [1 3]);
    
    %exit lightField
    exitOrigin = ppsfRays.origin;
    exitDirection = ppsfRays.direction;
    
    
    % do we need to put this in terms of angles? or are direction using
    % cartesian coordinates good enough?
    
    % put this into Ax = b form
    % A matrix will be a 4x4 matrix?   is it still 4x4 for 3d? maybe 5x5?
    % 6x6?
    % x will be a 4 x numSamples matrix containing the input lightfield
    % b will be a 4 x numSamples matrix containing the output lightfield
    % A will be the 4 x 4 least squares fit for this transformation
    
    % how to solve for A?  pseudo-inverse? SVD?
    
    
    
    
    
    
    
    %modify the rays for any aperture changes here
    modifyRays = ppsfObject();
    modifyRays.makeDeepCopy(ppsfRays);

    %trace from end of lens to sensor
    modifyRays.recordOnFilm(ppsfCamera.film);

    %show image
    ppsfCamera.showFilm();

    % save ppsf
    firstRays = ppsfObject();
    firstRays.makeDeepCopy(modifyRays);



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


