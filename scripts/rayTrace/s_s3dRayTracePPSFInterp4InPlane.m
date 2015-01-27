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
    
%% --first PPSF bottom right--
    %% point sources
    pointSources = lowerRightPosition;  %large distance test
    % pointSources = [ 0 0 -60];  %short distance test

    %% film properties
    film = filmC('position', [0 0 40 ],'size', [10 10], 'wave', 400:10:700);   %large distance

    %% lens properties
    %     diffractionEnabled = false;   
    lensFileName = fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat');
    lens = lensC('apertureSample', [301 301], 'fileName', lensFileName, 'apertureMiddleD', 16);
    lens.draw();

    %% ray trace and save ppsf
    ppsfCamera = ppsfCameraC('lens', lens, 'film', film, 'pointSource', pointSources);
    lowerRightRays = ppsfCamera.estimatePPSF();
    ppsfCamera.recordOnFilm();   % Could be on the camera, not the ppsf
    ppsfCamera.showFilm();
    
%% --2nd PPSF lower left --
    %% point sources
    pointSources = lowerLeftPosition;  %large distance test
    % pointSources = [ 0 0 -60];  %short distance test

    %% film properties -
    film = filmC('position', [0 0 40 ],'size', [10 10], 'wave', 400:10:700);   %large distance

    %% lens properties
    %     diffractionEnabled = false;   
    lensFileName = fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat');
    lens = lensC('apertureSample', [301 301], 'fileName', lensFileName, 'apertureMiddleD', 16);
    lens.draw();
    
    %% ray trace and save ppsf
    ppsfCamera = ppsfCameraC('lens', lens, 'film', film, 'pointSource', pointSources);
    lowerLeftRays = ppsfCamera.estimatePPSF();
    ppsfCamera.recordOnFilm();   % Could be on the camera, not the ppsf
    ppsfCamera.showFilm();
    
%% --3rd PPSF upper right--
    %% point sources
    pointSources = upperRightPosition;  %large distance test
    % pointSources = [ 0 0 -60];  %short distance test

    %% film properties -
    film = filmC('position', [0 0 40 ],'size', [10 10], 'wave', 400:10:700);   %large distance

    %% lens properties
    lensFileName = fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat');
    lens = lensC('apertureSample', [301 301], 'fileName', lensFileName, 'apertureMiddleD', 16);
    lens.draw();
    
    %% ray trace and save ppsf
    ppsfCamera = ppsfCameraC('lens', lens, 'film', film, 'pointSource', pointSources);
    upperRightRays = ppsfCamera.estimatePPSF();
    ppsfCamera.recordOnFilm();   % Could be on the camera, not the ppsf
    ppsfCamera.showFilm();
    
%% --4th PPSF upper left--
    %% point sources
    pointSources = upperLeftPosition;  %large distance test
 
    %% film properties -
    film = filmC('position', [0 0 40 ],'size', [10 10], 'wave', 400:10:700);   %large distance

    %% lens properties
    lensFileName = fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat');
    lens = lensC('apertureSample', [301 301], 'fileName', lensFileName, 'apertureMiddleD', 16);
    lens.draw();
    
    %% ray trace and save ppsf
    ppsfCamera = ppsfCameraC('lens', lens, 'film', film, 'pointSource', pointSources);
    upperLeftRays = ppsfCamera.estimatePPSF();
    ppsfCamera.recordOnFilm();   % Could be on the camera, not the ppsf
    ppsfCamera.showFilm();   
    
%% --5th PPSF middle--
    %% point sources
    pointSources = newPosition;  %large distance test

       %% film properties -
    film = filmC('position', [0 0 40 ],'size', [10 10], 'wave', 400:10:700);   %large distance

    %% lens properties
    lensFileName = fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat');
    lens = lensC('apertureSample', [301 301], 'fileName', lensFileName, 'apertureMiddleD', 16);
    lens.draw();
    
    %% ray trace and save ppsf
    ppsfCamera = ppsfCameraC('lens', lens, 'film', film, 'pointSource', pointSources);
    middleRays = ppsfCamera.estimatePPSF();
    ppsfCamera.recordOnFilm();   % Could be on the camera, not the ppsf
    ppsfCamera.showFilm();
    
%% --Interpolate between 4 PPSFs to try to produce 5th one--
    %% Average between 4 sets of rays
    
    %perhaps make this into a function for ppsfrays....
    %ppsfRays.interpolate(point1, ppsf1, point2, ppsf2, point3, ppsf3, point4, ppsf4);
    
    averageRays = ppsfC();
    averageRays.makeDeepCopy(lowerRightRays);
    %average between the first 2 ppsfs
    averageRays.origin = (lowerRightWeight * lowerRightRays.origin + lowerLeftWeight* lowerLeftRays.origin + upperLeftWeight * upperLeftRays.origin + upperRightWeight * upperRightRays.origin);
    averageRays.direction = (lowerRightWeight * lowerRightRays.direction + lowerLeftWeight* lowerLeftRays.direction + upperLeftWeight* upperLeftRays.direction + upperRightWeight*upperRightRays.direction);

    
    %other properties that must be interpolated...
    averageRays.aEntranceInt;
    averageRays.aMiddleInt;
    averageRays.aExitInt;
    averageRays.aExitDir;
    
    
    
    %     averageRays.apertureLocation = (lowerRightWeight * lowerRightRays.apertureLocation + lowerLeftWeight * lowerLeftRays.apertureLocation + upperLeftWeight * upperLeftRays.apertureLocation + upperRightWeight* upperRightRays.apertureLocation);

    % ray-trace the last bit - from lens to sensor
    %% modify the film and see the consequences on the PSF - these computations
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
%     modifyRays.origin(outsideAperture, : ) = [];   %this needs to be fixed later
%     modifyRays.direction(outsideAperture, : ) = [];
%     modifyRays.wavelength(outsideAperture) = [];
%     modifyRays.waveIndex(outsideAperture) = [];
%     modifyRays.apertureLocation(outsideAperture, :) = [];
%     modifyRays.apertureSamples.X(outsideAperture) = []; 
%     modifyRays.apertureSamples.Y(outsideAperture) = [];
    
    %% record on film
    filmCell = cell(1,1);
    %first try at 36.4 sensor distance
    film = filmC('position', [0 0 40 ],'size', [10 10], 'wave', 400:10:700);   %large distance
    filmCell{1} = film;
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


