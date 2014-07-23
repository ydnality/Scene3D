


%% use only 1 point source for now - raytrace through the lens

%calculate the origin and direction of the rays
%     rays.traceSourceToLens(pointSources(curInd, :), lens);

disp('-----trace source to lens-----');
tic
rays = lens.rayTraceSourceToLens(pointSources(1, :));

apertureSamples = lens.apertureSample;
ppsfRays = ppsfC(rays.origin, rays.direction, rays.wavelength, pointSourceDepth, pointSourceFieldHeight, apertureSamples);  %think of a best way to put in aperture sample location
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
modifyRays = ppsfC();
modifyRays.makeDeepCopy(ppsfRays);
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
filmCell{1} = film; %pbrtFilmC([0 0 36.4],[1 1], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance
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
