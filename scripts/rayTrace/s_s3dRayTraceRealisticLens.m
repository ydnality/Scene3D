%% ray-tracing for realistic lens
%
%  This uses Snell's law and a lens prescription to create a tray trace.
%  The script is too long, and we need to start writing functions so that
%  the length is shortened and the clarity increased.
%  We are only ray-tracing ideal point sources in order to extract out point
%  spread functions.
%
% AL Vistalab, 2014

%% ray-tracing 

% -New support: different wavelength support

% declare point sources in world space.  The camera is usually at [0 0 0],
% and pointing towards -z.  We are using a right-handed coordinate system.

%Large Distance Test
%place point sources using meshgrid.  This is a 5x5 grid.
% [XGrid YGrid] = meshgrid(-4000:1000:4000,-4000:1000:4000);   %large distance 
% pointSources = [XGrid(:) YGrid(:) ones(size(XGrid(:))) * -20000];
%only 1 point source in the middle
% pointSources = [ 0 0 -20000];  %large distance 

%Small Distance Test
[XGrid YGrid] = meshgrid(-15:5:15,-15:5:15);   %small distance test
pointSources = [XGrid(:) YGrid(:) ones(size(XGrid(:))) * -100];   %small distance ----ADJUST ME!!----
% pointSources = [ 0 0 -100];  %small distance - TURN THIS ON FOR ONLY 1

wave = [400 550 700];  % in nm
wavelengthConversion = [400 3; 550 2; 700 1];

%% film properties - 
film = filmObject([], [],  400:10:700, [(400:10:700)' (1:31)'], []);

%% Should be a function for reading and writing lens files

% declare lens that has 2 simple elements
offset = [3, 0];
radius = [-67 67];
aperture = [3 3];
n = [ 1 1.67];
lensCenterPosition = [0 0 -1.5];  %eventually calculate this given the lens file
lens = lensRealisticObject(offset,radius,aperture,n, 3, lensCenterPosition);

%% loop through all point sources
vcNewGraphWin; %for illustration
for curInd = 1:size(pointSources, 1);
    %calculate the origin and direction of the rays
    rays = rayObject;
    %     rays.traceSourceToLens(pointSources(curInd, :), lens);
    lens.rayTraceSourceToLens(pointSources(curInd, :), rays);
    
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
vcAddAndSelectObject(oi); oiWindow;

%% End
