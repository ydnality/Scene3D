%% ray-tracing for realistic lens
%
%  This uses Snell's law and a lens prescription to create a tray trace.
%  The script is too long, and we need to start writing functions so that
%  the length is shortened and the clarity increased.
%  We are only ray-tracing ideal point sources in order to extract out point
%  spread functions.
%
% AL Vistalab, 2014
%%
s_initISET
%% point sources
pointSources = [ 0 0 -20000];  %large distance test
% pointSources = [ 0 0 -60];  %short distance test

%% film properties -
film = filmC([0 0 36.4],[1 1], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance

%% lens properties
diffractionEnabled = true;

%initialize to default
%lensRealisticObject(elOffset, elRadius, elAperture, elN, aperture, focalLength, center, diffractionEnabled, wave)
lens = lensRealisticObject([],[],[],[], 8, 50, [], diffractionEnabled);

%read lens from file
lens.readLensFile(fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat'))

%intialize ray=tracing
lens.calculateApertureSample([201 201]);

%lens illustration
lens.drawLens();

%% loop through all point sources
%TODO - put this into function - perhaps a "camera" object?
% vcNewGraphWin; %for illustration
for curInd = 1:size(pointSources, 1);
    %calculate the origin and direction of the rays
    %     rays.traceSourceToLens(pointSources(curInd, :), lens);
    
    disp('-----trace source to lens-----');
    tic
    rays = lens.rayTraceSourceToLens(pointSources(curInd, :));
    toc
    
    %duplicate the existing rays, and creates one for each
    %wavelength
    disp('-----expand wavelenghts-----');
    tic
    rays.expandWavelengths(film.wave);
    toc
    
    %lens intersection and raytrace
    disp('-----rays trace through lens-----');
    tic
    lens.rayTraceThroughLens(rays);
    toc
    
    %intersect with "film" and add to film
    disp('-----record on film-----');
    tic
    rays.recordOnFilm(film);
    toc
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

