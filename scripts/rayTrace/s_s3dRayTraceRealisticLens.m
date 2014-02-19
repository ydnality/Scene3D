%% ray-tracing for realistic lens
%
%  This uses Snell's law and a lens prescription to create a tray trace.
%  The script is too long, and we need to start writing functions so that
%  the length is shortened and the clarity increased.
%
% AL Vistalab, 2014

%% ray-tracing 
% -We are only ray-tracing ideal point sources in order to extract out point
% spread functions.
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

%% sensor properties - 
% sensor.position = [0 0 49];  %large distance 
% sensor.position = [0 0 165]; 
sensor.position = [0 0 100];  %small distance ----ADJUST ME!!----
sensor.size = [48 48];  %in mm
sensor.resolution = [200 200 length(wave)];
sensor.image = zeros(sensor.resolution);

% sensor.focalLength = 50; %in mm
% apertureRadius = 1; % in mm
apertureRadius =3; % in mm


%% Should be a function for reading and writing lens files

%note: the offset format is DIFFERENT from the lens files - we must fix
%this somehow (AL)
lensSurfaces = cell(1,1);
lensSurfaces{1}.offset = 3;
lensSurfaces{1}.radius = -67;
lensSurfaces{1}.aperture = 3;
lensSurfaces{1}.indexOfRefr = 1;

lensSurfaces{2}.offset = 0;   %figure out the best way to deal with these offsets
lensSurfaces{2}.radius = 67;
lensSurfaces{2}.aperture = 3;
lensSurfaces{2}.indexOfRefr = 1.67;

% Start testing here ...
%
% offset = [3, 0]
% radius = [-67 67]
% aperture = [3 3]
% n = [ 1 1.67];
% lensObject(offset,radius,aperture,n)

totalOffset = 0;
for i = 1:length(lensSurfaces)
    totalOffset = totalOffset + lensSurfaces{i}.offset; 
end

%% Set up aperture sampling positions 

% This should be a function
% Create aperture sample positions
% Adjust the sampling rate if needed - this determines the number of
% samples per light source 

%loop through aperture positions and uniformly sample the aperture
%everything is done in vector form for speed
[apertureSample.X, apertureSample.Y] = meshgrid(linspace(-1, 1, 20),linspace(-1, 1, 20)); %adjust this if needed - this determines the number of samples per light source

% this position should NOT change
lensCenterPosition = [0 0 -1.5];

%debug illustrations initialize
lensIllustration = zeros(300, 300);

%% loop through all point sources
vcNewGraphWin;

for curInd = 1:size(pointSources, 1);
    curPointSource = pointSources(curInd, :);

    %assume a circular aperture, and make a mask that is 1 when the pixel
    %is within a circle of radius 1
    apertureMask = (apertureSample.X.^2 + apertureSample.Y.^2) <= 1;
    
    scaledApertureSample.X = apertureSample.X .* apertureRadius; 
    scaledApertureSample.Y = apertureSample.Y .* apertureRadius; 
    
    %remove cropped sections of aperture
    croppedApertureSample.X =  scaledApertureSample.X;
    croppedApertureSample.X(apertureMask == 0) = [];
    croppedApertureSample.X = croppedApertureSample.X';
    croppedApertureSample.Y =  scaledApertureSample.Y;
    croppedApertureSample.Y(apertureMask == 0) = [];
    croppedApertureSample.Y = croppedApertureSample.Y';
    
    %calculate the origin and direction of the rays
    rays.origin = repmat(curPointSource, [size(croppedApertureSample.Y(:), 1) 1] );   %the new origin will just be the position of the current light source
    rays.direction = [(croppedApertureSample.X(:) -  rays.origin(:,1)) (croppedApertureSample.Y(:) -  rays.origin(:,2)) (lensCenterPosition(3) - rays.origin (:,3)) .* ones(size(croppedApertureSample.Y(:)))];
    rays.direction = rays.direction./repmat( sqrt(rays.direction(:, 1).^2 + rays.direction(:, 2).^2 + rays.direction(:,3).^2), [1 3]); %normalize direction

    %add different wavelengths
    %first duplicate the existing entries, and create one for each
    %wavelength
    subLength = size(rays.origin, 1);
    rays.origin = repmat(rays.origin, [length(wave) 1]);
    rays.direction = repmat(rays.direction, [length(wave) 1]);
    
    %creates a vector representing wavelengths... for example: [400 400 400... 410 410 410... ..... 700]
    tmp = (wave' * ones(1, subLength))'; 
    rays.wavelength = tmp(:);  
    
    
    prevSurfaceZ = -totalOffset;
    prevN = 1;  %assume that we start off in air

    %initialize newRays to be the old ray.  We will update it later.
    newRays = rays;
    for lensEl = length(lensSurfaces):-1:1
        sphere = lensSurfaces{lensEl};
        sphere.center = [0 0 prevSurfaceZ + sphere.offset + sphere.radius];
        
        %illustrations for debug
        zPlot = linspace(sphere.center(3) - sphere.radius, sphere.center(3) + sphere.radius, 10000);
        yPlot = sqrt(sphere.radius^2 - (zPlot - sphere.center(3)) .^2);
        yPlotN = -sqrt(sphere.radius^2 - (zPlot - sphere.center(3)) .^2);
        arcZone = 5;
        %TODO:find a better way to plot the arcs later - this one is prone to potential problem
        withinRange = and(and((yPlot < sphere.aperture),(zPlot < prevSurfaceZ + sphere.offset + arcZone)), (zPlot > prevSurfaceZ + sphere.offset - arcZone));  
        line(zPlot(withinRange), yPlot(withinRange));
        line(zPlot(withinRange), yPlotN(withinRange));   
        
        %vectorize this operation later
        for i = 1:size(rays.origin, 1)
            %get the current ray
            ray.direction = newRays.direction(i,:);
            ray.origin = newRays.origin(i,:);
            ray.wavelength = newRays.wavelength(i);
            
            %calculate intersection with spherical lens element
            radicand = dot(ray.direction, ray.origin - sphere.center)^2 - ...
                ( dot(ray.origin -sphere.center, ray.origin -sphere.center)) + sphere.radius^2;
            if (sphere.radius < 0)
                intersectT = (-dot(ray.direction, ray.origin - sphere.center) + sqrt(radicand));
            else
                intersectT = (-dot(ray.direction, ray.origin - sphere.center) - sqrt(radicand));
            end

            %make sure that T is > 0
            if (intersectT < 0)
                disp('Warning: intersectT less than 0.  Something went wrong here...');
            end
            
            intersectPosition = ray.origin + intersectT * ray.direction;
            
            %illustrations for debugging
            lensIllustration(max(round(intersectPosition(2) * 100 + 150),1), max(-round(intersectPosition(3) * 1000), 1)) = 1;  %show a lens illustration
            hold on;
            line([ray.origin(3) intersectPosition(3) ], [ray.origin(2) intersectPosition(2)] ,'Color','b','LineWidth',1);
            

            normalVec = intersectPosition - sphere.center;  %does the polarity of this vector matter? YES
            normalVec = normalVec./norm(normalVec);
            if (sphere.radius < 0)  %which is the correct sign convention? This is correct
                normalVec = -normalVec;
            end
            
            
            %modify the index of refraction depending on wavelength
            %TODO: have this be one of the input parameters (N vs. wavelength)
            if (sphere.indexOfRefr ~= 1)
                curN = (ray.wavelength - 550) * -.04/(300) + sphere.indexOfRefr;
            else
                curN = 1;
            end
            
            ratio = prevN/curN;    %snell's law index of refraction

            %Vector form of Snell's Law
            c = -dot(normalVec, ray.direction);
            newVec = ratio *ray.direction + (ratio*c -sqrt(1 - ratio^2 * (1 - c^2)))  * normalVec; 
            newVec = newVec./norm(newVec); %normalize
            
            %update the direction of the ray
            newRays.origin(i, : ) = intersectPosition;
            newRays.direction(i, : ) = newVec;
        end
        prevN = curN;
        
        
        prevSurfaceZ = prevSurfaceZ + sphere.offset;
     end

    %calculate intersection point at sensor
    intersectZ = repmat(sensor.position(3), [size(newRays.origin, 1) 1]);
    intersectT = (intersectZ - newRays.origin(:, 3))./newRays.direction(:, 3);
    intersectPosition = newRays.origin + newRays.direction .* repmat(intersectT, [1 3]);
    
    
    %imagePixel is the pixel that will gain a photon due to the traced ray
    imagePixel.position = [intersectPosition(:,2) intersectPosition(:, 1)]; 
    imagePixel.position = real(imagePixel.position); %add error handling for this
    imagePixel.position = round(imagePixel.position * sensor.resolution(1)/sensor.size(1) + ...
        repmat( sensor.resolution(1:2)./2, [size(imagePixel.position,1) 1]));   %
    %scale the position to a sensor position
    imagePixel.position(imagePixel.position < 1) = 1; %make sure pixel is in range
    imagePixel.position = min(imagePixel.position, repmat(sensor.resolution(1:2), [size(imagePixel.position,1) 1]));
    imagePixel.wavelength = newRays.wavelength; 
    
    
    
    %add a value to the intersection position
    for i = 1:size(rays.origin , 1)
        wantedPixel = [imagePixel.position(i,1) imagePixel.position(i,2) find(wavelengthConversion == imagePixel.wavelength(i))];  %pixel to update
        sensor.image(wantedPixel(1), wantedPixel(2), wantedPixel(3)) =  sensor.image(wantedPixel(1), wantedPixel(2), wantedPixel(3)) + 1;  %sensor.image(imagePixel(:,1), imagePixel(:,2)) + 1;
        %illustrations for debugging
        line([newRays.origin(i, 3) intersectPosition(i, 3)] ,  [newRays.origin(i, 2);  intersectPosition(i, 2)]);
    end
end

%% Show the image

vcNewGraphWin;
imshow(sensor.image/ max(sensor.image(:)));

%% End
