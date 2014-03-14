fid = fopen(fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat'));
import = textscan(fid, '%s%s%s%s', 'delimiter' , '\t');
fclose(fid);


%first find the start of the lens line demarked by #   radius

firstColumn = import{1};
startLens = false;
elementCount = 1;


%read comment lines
continu = true;
i = 1;
while(continu && i <= length(firstColumn))
    compare = regexp(firstColumn(i), 'radius');
    if(~(isempty(compare{1})))
      continu = false; 
    end
    i = i+1;
end

%put data into lens object
radius = str2double(import{1});  
radius = radius(i:length(firstColumn));

offset = str2double(import{2});
offset = offset(i:length(firstColumn));
%change from pbrt Scene3D format to raytrace Scene3D format
offset = [0; offset(1:(end-1))];


N = str2double(import{3});
N = N(i:length(firstColumn));

aperture = str2double(import{4})/2;%radius supplied is the radius diameter, so divide it by 2
aperture = aperture(i:length(firstColumn));


%% end read file

pointSources = [ 0 0 -20000]; 
film = filmObject([0 0 36.4],[1 1], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance



lensCenterPosition = [0 0 -15.1550];  %eventually calculate this given the lens file

diffractionEnabled = false;
% lensRealisticObject(elOffset, elRadius, elAperture, elN, aperture, focalLength, center, diffractionEnabled)
lens = lensRealisticObject(offset,radius,aperture,N, 1, 50, lensCenterPosition, diffractionEnabled);
lens.calculateApertureSample([51 51]);
lens.drawLens();

%% loop through all point sources
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

