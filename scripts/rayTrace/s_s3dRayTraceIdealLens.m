%% Ray-tracing for an ideal lens
%
% An ideal thin lens is defined as one that obey's the thin lens equation
% (1/s1 + 1/s2 = 1/f).  Therefore, Snell's Law is not used in this script
% at all.  Instead, the direction that the rays bend are determined by the
% thin lens equation.  At each point source, a "center ray" is shot at the
% center of the lens.  The intersection of this ray and the focal-plane as
% defined by the thin lens equation determines the point of focus of the
% lens.  All other rays that are shot at the edge of the aperture will then
% intersect this ideal focal point.  
% All units are in mm.
%
%  TODO:
%   We can't use sensor here and sensorCreate in same way.
%   Somehow we should integrate this sensor stuff with ISET.
%   I don't see the lens part here.  Need more comments about what this is
%   doing, and what it tests.
%
% See also: s_3dRayTrace*.m
%
% AL Vistalab 2014

diffractionEnabled = true;
%% point sources
% declare point sources in world space.  The camera is usually at [0 0 0],
% and pointing towards -z.  We are using a right-handed coordinate system.

% [XGrid YGrid] = meshgrid(-4000:1000:4000,-4000:1000:4000);
% pointSources = [XGrid(:) YGrid(:) ones(size(XGrid(:))) * -20000];
pointSources = [0 0 -20000];
 


%% camera properties 
% TODO: change the name
% sensor.size = [1 1]; 
% %sensor.size = [48 48]; %in mm
% sensor.position = [0 0 50]; 
% sensor.resolution = [200 200 length(wave)];
% sensor.image = zeros(sensor.resolution);
% sensor.focalLength = 50; %in mm
% apertureRadius =.1; % in mm



%Create a set of circular aperture positions and uniformly sample the aperture
%everything.  Calculations will be done in vector form for speed
% The code here is a good candidate for a function (BW).
% [apertureSample.X, apertureSample.Y] = meshgrid(linspace(-1, 1, 90),linspace(-1, 1, 90)); %adjust this if needed


% this position should NOT change
% lensCenterPosition = [0 0 0];
film = filmObject([0 0 50],[1 1], 400:10:700, [(400:10:700)' (1:31)'], []);   %(position, size,  wave, waveConversion, resolution)
lens = lensIdealObject(.1, 50, [0 0 0], true);

lens.calculateApertureSample([20 20]);


%% loop through all point sources
for curInd = 1:size(pointSources, 1);
    
    % This calculation happens a lot ... we should functionalize it.
%      curPointSource = pointSources(curInd, :);
    
    

    %calculate the origin and direction of the rays
%     rays.origin = repmat(curPointSource, [size(lens.apertureSample.Y(:), 1) 1] );   %the new origin will just be the position of the current light source
%     rays.direction = [(lens.apertureSample.X(:) -  rays.origin(:,1)) (lens.apertureSample.Y(:) -  rays.origin(:,2)) (lens.centerPosition(3) - rays.origin (:,3)) .* ones(size(lens.apertureSample.Y(:)))];
%     rays.direction = rays.direction./repmat( sqrt(rays.direction(:, 1).^2 + rays.direction(:, 2).^2 + rays.direction(:,3).^2), [1 3]); %normalize direction
%        
%     
    rays = rayObject;
    lens.rayTraceSourceToLens(pointSources(curInd, :), rays);
    

    %duplicate the existing rays, and creates one for each
    %wavelength
    rays.expandWavelengths(film.wave);
    
%     %first duplicate the existing entries, and create one for each
%     %wavelength
%     subLength = size(rays.origin, 1);
%     rays.origin = repmat(rays.origin, [length(film.wave) 1]);
%     rays.direction = repmat(rays.direction, [length(film.wave) 1]);
%     %creates a vector representing wavelengths... for example: [400 400 400... 410 410 410... ..... 700]
%     rays.wavelength = vectorize((film.wave' * ones(1, subLength))');  
%     
%     
    
%--lens refraction
    
    lens.rayTraceThroughLens(rays, pointSources(curInd, :));


%     
% %     %calculate intersection point at sensor
%     intersectZ = repmat(film.position(3), [size(newRays.origin, 1) 1]);
%     intersectT = (intersectZ - newRays.origin(:, 3))./newRays.direction(:, 3);
%     intersectPosition = newRays.origin + newRays.direction .* repmat(intersectT, [1 3]);
%     
%      %imagePixel is the pixel that will gain a photon due to the traced ray
%     imagePixel.position = [intersectPosition(:,2) intersectPosition(:, 1)]; 
%     imagePixel.position = real(imagePixel.position); %add error handling for this
%     imagePixel.position = round(imagePixel.position * film.resolution(1)/film.size(1) + ...
%         repmat( film.resolution(1:2)./2, [size(imagePixel.position,1) 1]));   %
%     %scale the position to a film position
%     imagePixel.position(imagePixel.position < 1) = 1; %make sure pixel is in range
%     imagePixel.position = min(imagePixel.position, repmat(film.resolution(1:2), [size(imagePixel.position,1) 1]));
%     imagePixel.wavelength = newRays.wavelength; 
%     
% 
%     %add a value to the intersection position
%     for i = 1:size(rays.origin , 1)
%         wantedPixel = [imagePixel.position(i,1) imagePixel.position(i,2) find(film.waveConversion == imagePixel.wavelength(i))];  %pixel to update
%         film.image(wantedPixel(1), wantedPixel(2), wantedPixel(3)) =  film.image(wantedPixel(1), wantedPixel(2), wantedPixel(3)) + 1;  %sensor.image(imagePixel(:,1), imagePixel(:,2)) + 1;
%         %illustrations for debugging
% %         line([newRays.origin(i, 3) intersectPosition(i, 3)] ,  [newRays.origin(i, 2);  intersectPosition(i, 2)]);
%     end


%     intersect with "film" and add to film
   rays.recordOnFilm(film);

end

%%
% vcNewGraphWin; imshow(sensor.image/ max(sensor.image(:)));


%assign as optical image
oi = oiCreate;
oi = initDefaultSpectrum(oi);
oi = oiSet(oi,'photons',film.image);
vcAddAndSelectObject(oi); oiWindow;


%%
