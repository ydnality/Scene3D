%% ray-tracing for realistic lens - PPSF
%
% Lens element positions are all negative, and the final exit plan can be
% considered as z = 0
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
%
% AL Vistalab, 2014
%%
s_initISET

%% --first PPSF at the center--
    %% point sources
    pointSourceDepth = 1000;
    pointSources = [ 0 0 -pointSourceDepth];  %large distance test
    pointSourceFieldHeight = 0;
    % pointSources = [ 0 0 -60];  %short distance test

    %% film properties
    % filmPosition = [ 0 0 40];
    % filmSize = [10 10]; % mm
    % wave = 400:50:700;
    % next will go
    % don't remember
    
    %40
    film = pbrtFilmObject([0 0 60 ],[10 10], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance

    %% lens properties
    diffractionEnabled = false;   
    %turning on diffraction does NOT make sense yet since we have not modeled 
    %the transformation of uncertainty from the middle aperture to the end of the lens

    %initialize to default
    %lensRealisticObject(elOffset, elRadius, elAperture, elN, aperture, focalLength, center, diffractionEnabled, wave)
    lens = lensMEObject('apertureSample', [201 201]);

    %read lens from file
    lens.fileRead(fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat'))

    %lens illustration
    lens.drawLens();

    %% ray trace and save ppsf - Not sure camera should have pointSources
    ppsfCamera = ppsfCameraObject('lens', lens, 'film', film, 'pointSource', pointSources);
    ppsf = ppsfCamera.estimatePPSF();
    
    %     %project the rays from the final lens onto the exit pupil plane at the
    %     %z = 0 (exit pupil plane).
    ppsf.projectOnPlane(0);
    
    %% Calculate light field at the entrance pupil plane and exit pupil - estimate the linear transform
    
%     % These are the X,Y samples in the entrance pupil, which corresponds to
%     % the most negative point(first surface) of the lens elements
%     initialOrigin = ppsf.aEntranceInt;   % Only has X and Y
%     
%     % The vector that connects the point on the pupil plane to the point
%     % source
%     initialDirection = zeros(size(initialOrigin.X), 3);
%     initialDirection(:,1) = initialOrigin.X - pointSources(1);
%     initialDirection(:,2) = initialOrigin.Y - pointSources(2);
%     initialDirection(:,3) = -lens.get('totaloffset') - pointSources(3);
%     
%     % Make this a make unit length function
%     initialDirection = initialDirection./repmat(sum(initialDirection,2), [1 3]);
%     
%     %exit lightField
%     exitOrigin = ppsf.aExitInt;
%     exitDirection = ppsf.aExitDir;
%     
%     
%     % do we need to put this in terms of angles? or are direction using
%     % cartesian coordinates good enough?
%     
%     % put this into Ax = b form
%     % A matrix will be a 4x4 matrix?   is it still 4x4 for 3d? maybe 5x5?
%     % 6x6?
%     % x will be a 4 x numSamples matrix containing the input lightfield
%     % b will be a 4 x numSamples matrix containing the output lightfield
%     % A will be the 4 x 4 least squares fit for this transformation
%     
%     % how to solve for A?  pseudo-inverse? SVD?
    
    
    %% record on film
    

    %modify the rays for any aperture changes here
    modifyRays = ppsfObject();
    modifyRays.makeDeepCopy(ppsfCamera.ppsfRays);

    %trace from end of lens to sensor
    modifyRays.recordOnFilm(ppsfCamera.film);
    %show image
    ppsfCamera.showFilm();

    % save ppsf
%     firstRays = ppsfObject();
%     firstRays.makeDeepCopy(modifyRays);


    %% Show the images
    
    % vcNewGraphWin;
    % imshow(film.image/ max(film.image(:)));
   
%         oi = oiCreate;
%         oi = initDefaultSpectrum(oi);
%         oi = oiSet(oi, 'wave', film.wave);
%         oi = oiSet(oi,'photons',film.image);
% 
% 
%         optics = oiGet(oi,'optics');
%         optics = opticsSet(optics,'focal length',lens.focalLength/1000);
%         optics = opticsSet(optics,'fnumber', lens.focalLength/(2*1));
%         oi = oiSet(oi,'optics',optics);
%         hfov = rad2deg(2*atan2(film.size(1)/2,lens.focalLength));
%         oi = oiSet(oi,'hfov', hfov);
% 
%         temp = film.position;
%         filmDistance = temp(3);
%         oi = oiSet(oi, 'name', ['filmDistance: ' num2str(filmDistance)]);
%         vcAddAndSelectObject(oi); oiWindow;


