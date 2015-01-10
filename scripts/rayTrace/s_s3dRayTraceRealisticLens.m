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
ieInit

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


%% film properties - 
film = pbrtFilmC('wave', 400:10:700);  %small distance

%% Should be a function for reading and writing lens files

% declare lens that has 2 simple elements
% offset = [1.5 1.5 0];
% radius = [-67 0 67];
% aperture = [5 1 4];
% n = [ 1 0 1.67];

offset = [0 1.5 1.5];
radius = [67 0 -67];
aperture = [4 4 5];
% aperture = [5 .5 1];

n = [ 1.67 0 1];
lensCenterPosition = [0 0 -1.5];  %eventually calculate this given the lens file

diffractionEnabled = true;   %reenable
diffractionEnabled = false;

lens = lensC; % (offset,radius,aperture,n, 8, 50, lensCenterPosition, diffractionEnabled);
lens.calculateApertureSample([7 7]);

%% loop through all point sources
vcNewGraphWin; %for illustration
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

vcAddObject(oi); oiWindow;


%TODO: allow for real lens aperture, not just a sampling function change V
%(make rays terminate when they hit the black section of aperture

%TODO: vectorize rays  V
%TODO: see if we can create a more efficient aperture sampling procedure,
%rather than sampling from the furthest most aperture

%TODO: fix film recording action if it is out of bounds    v
%TODO: try to figure out what is the deal with the non-real numbers
%produced 


%% End
