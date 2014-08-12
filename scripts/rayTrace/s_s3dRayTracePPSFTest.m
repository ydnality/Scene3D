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
% AL Vistalab, 2014
%%
s_initISET
%% point sources
pointSourceDepth = 20000;
pointSourceFieldHeight = 0;
pointSources = [ 3 0 -pointSourceDepth];  %large distance test
% pointSources = [ 0 0 -60];  %short distance test

%% film properties -
film = pbrtFilmC([0 0 36.4],[1 1], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance

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
randJitter = true;
lens.calculateApertureSample([301 301], randJitter);

%lens illustration
lens.drawLens();

%% use only 1 point source for now - raytrace through the lens

    
%calculate the origin and direction of the rays
%     rays.traceSourceToLens(pointSources(curInd, :), lens);

disp('-----trace source to lens-----');
tic
rays = lens.rayTraceSourceToLens(pointSources(1, :));

apertureSamples = lens.apertureSample;
ppsfRays = ppsfC('origin', rays.origin, 'direction', rays.direction, ...
    'waveIndex', rays.waveIndex, 'pointSourceDepth', pointSourceDepth, ...
    'pointSourceFieldHeight', pointSourceFieldHeight, ...
    'aEntranceInt', apertureSamples);  %think of a best way to put in aperture sample location
toc

%duplicate the existing rays, and creates one for each
%wavelength
disp('-----expand wavelenghts-----');
tic
ppsfRays.expandWavelengths(film.wave);
toc

%lens intersection and raytrace
disp('-----rays trace through lens-----');
tic
lens.rayTraceThroughLens(ppsfRays);
toc

%% The rays at this point can then be saved to file and stored as precomputed rays

%TODO: ray saving, and ray loading

%% ray-trace the last bit - from lens to sensor
%modify the film and see the consequences on the PSF - these computations
%should be very fast
modifyRays = ppsfObject();
modifyRays.makeDeepCopy(ppsfRays);
% 
newRadius = 2;
outsideAperture = modifyRays.apertureLocation(:,1).^2 + modifyRays.apertureLocation(:,2).^2 > newRadius^2;

%modify so only x >0 shows up
% outsideAperture = modifyRays.apertureSamples.X > 0; 

%remove outside of aperture elements
%TODO: make this into a function
modifyRays.origin(outsideAperture, : ) = [];
modifyRays.direction(outsideAperture, : ) = [];
%modifyRays.wavelength(outsideAperture) = [];
modifyRays.waveIndex(outsideAperture) = [];
modifyRays.apertureLocation(outsideAperture, :) = [];
modifyRays.apertureSamples.X(outsideAperture) = []; 
modifyRays.apertureSamples.Y(outsideAperture) = [];

film = cell(1,3);
%first try at 36.4 sensor distance
film{1} = pbrtFilmC([0 0 36.4],[1 1], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance
%intersect with "film" and add to film
disp('-----record on film-----');
tic
modifyRays.recordOnFilm(film{1});
toc

%38 sensor distance
film{2} = pbrtFilmC([0 0 37],[1 1], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance
%intersect with "film" and add to film
disp('-----record on film-----');
tic
modifyRays.recordOnFilm(film{2});
toc

%35 sensor distance
film{3} = pbrtFilmC([0 0 35.5],[1 1], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance
%intersect with "film" and add to film
disp('-----record on film-----');
tic
modifyRays.recordOnFilm(film{3});
toc

%% Show the images

% vcNewGraphWin;
% imshow(film.image/ max(film.image(:)));

for i = 1:length(film)

    oi = oiCreate;
    oi = initDefaultSpectrum(oi);
    oi = oiSet(oi, 'wave', film{i}.wave);
    oi = oiSet(oi,'photons',film{i}.image);


    optics = oiGet(oi,'optics');
    optics = opticsSet(optics,'focal length',lens.focalLength/1000);
    optics = opticsSet(optics,'fnumber', lens.focalLength/(2*1));
    oi = oiSet(oi,'optics',optics);
    hfov = rad2deg(2*atan2(film{i}.size(1)/2,lens.focalLength));
    oi = oiSet(oi,'hfov', hfov);
    
    temp = film{i}.position;
    filmDistance = temp(3);
    oi = oiSet(oi, 'name', ['filmDistance: ' num2str(filmDistance)]);
    vcAddAndSelectObject(oi); oiWindow;
end