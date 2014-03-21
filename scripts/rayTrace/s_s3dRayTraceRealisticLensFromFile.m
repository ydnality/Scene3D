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

%% film properties - 
film = filmObject([0 0 36.4],[1 1], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance

%% lens properties
diffractionEnabled = false;

%initialize to default
lens = lensRealisticObject([],[],[],[], 2, 50, [], diffractionEnabled);

%read lens from file
lens.readLensFile(fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat'))

%intialize ray=tracing
lens.calculateApertureSample([101 101]);

%lens illustration
lens.drawLens();

%% loop through all point sources
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

