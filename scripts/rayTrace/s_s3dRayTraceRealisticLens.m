%% ray-tracing for an ideal lens

% declare point sources in world space.  The camera is usually at [0 0 0],
% and pointing towards -z.  We are using a right-handed coordinate system.

%place point sources using meshgrid.  This is a 5x5 grid.
[XGrid YGrid] = meshgrid(-4000:1000:4000,-4000:1000:4000);
pointSources = [XGrid(:) YGrid(:) ones(size(XGrid(:))) * -20000];

% [XGrid YGrid] = meshgrid(-400:100:400,-400:100:400);
% pointSources = [XGrid(:) YGrid(:) ones(size(XGrid(:))) * -7000];

% pointSources = [0 0 -20000;
%                 0 5000 -20000;
%                 0 -5000 -20000;
%                 5000 0 -20000;
%                 5000 5000 -20000;
%                 5000 -5000 -20000;
%                 -5000 0 -20000;
%                 -5000 5000 -20000;
%                 -5000 -5000 -20000;];
 
            
            
            

%sensor properties
sensor.size = [48 48];  %in mm
sensor.position = [0 0 49]; 
% sensor.position = [0 0 165]; 
sensor.resolution = [200 200];
sensor.image = zeros(sensor.resolution);
sensor.focalLength = 50; %in mm

apertureRadius =1; % in mm

%note: the offset format is DIFFERENT from the lens files - we must address
%this somehow
lensSurfaces = cell(1,1);
lensSurfaces{1}.offset = 3;
lensSurfaces{1}.radius = -67;
lensSurfaces{1}.aperture = 3;
lensSurfaces{1}.indexOfRefr = 1;

lensSurfaces{2}.offset = 0;   %figure out the best way to deal with these offsets
lensSurfaces{2}.radius = 67;
lensSurfaces{2}.aperture = 3;
lensSurfaces{2}.indexOfRefr = 1.67;

totalOffset = 0;
for i = 1:length(lensSurfaces)
    totalOffset = totalOffset + lensSurfaces{i}.offset; 
end



% this position should NOT change
lensCenterPosition = [0 0 1.5];

lensIllustration = zeros(300, 300);
%loop through all point sources
for curInd = 1:size(pointSources, 1);
    
    curPointSource = pointSources(curInd, :);
    
    %trace ray from point source to lens center, to image.  This helps
    %determine the point of focus
    centerRay.origin = curPointSource;
    centerRay.direction = lensCenterPosition - centerRay.origin;
    centerRay.direction = centerRay.direction./norm(centerRay.direction);
    
    %calculate the in-focus position using thin lens equation
    inFocusDistance = 1/(1/sensor.focalLength - -1/curPointSource(3));
    inFocusT = (inFocusDistance - centerRay.origin(3))/centerRay.direction(3);
    inFocusPosition = centerRay.origin + inFocusT .* centerRay.direction;
    

    %loop through aperture positions and uniformly sample the aperture
    %everything is done in vector form for speed
    [apertureSample.X, apertureSample.Y] = meshgrid(linspace(-1, 1, 30),linspace(-1, 1, 30)); %adjust this if needed
    
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
       
    
    
    prevSurfaceZ = -totalOffset;
    prevN = 1;  %assume that we start off in air
    
    
    for lensEl = length(lensSurfaces):-1:1
        sphere = lensSurfaces{lensEl};
        sphere.center = [0 0 prevSurfaceZ + sphere.offset + sphere.radius];

        newRays = rays;
        %vectorize this operation later
        for i = 1:length(rays.origin)
            ray.direction = newRays.direction(i,:);
            ray.origin = newRays.origin(i,:);
            %calculate intersection with lens element
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
            lensIllustration(max(round(intersectPosition(2) * 100 + 150),1), max(-round(intersectPosition(3) * 1000), 1)) = 1;  %show a lens illustration
            
            normalVec = intersectPosition - sphere.center;  %does the polarity of this vector matter?
            normalVec = normalVec./norm(normalVec);
            if (sphere.radius < 0)  %which is the correct sign convention?
                normalVec = -normalVec;
            end
            ratio = prevN/sphere.indexOfRefr;
            c = -dot(normalVec, ray.direction);
            newVec = ratio *ray.direction + (ratio*c -sqrt(1 - ratio^2 * (1 - c^2)))  * normalVec; 
            newVec = newVec./norm(newVec); %normalize
            newRays.origin(i, : ) = intersectPosition;
            newRays.direction(i, : ) = newVec;
            
        end
        prevN = sphere.indexOfRefr;
        prevSurfaceZ = prevSurfaceZ + sphere.offset;
     end
%     %when intersecting ideal lens, change the direction to intersect the
%     %inFocusPosition, and update the origin
%     lensIntersectT = (lensCenterPosition(3) - rays.origin(:,3))./ rays.direction(:,3);
%     lensIntersectPosition = rays.origin +  repmat(lensIntersectT, [1 3]) .* rays.direction;
%     %calculate new direction
%     newRays.origin = lensIntersectPosition;
%     newRays.direction = repmat(inFocusPosition , [size(newRays.origin,1) 1 ]) - newRays.origin;
%     

    %calculate intersection point at sensor
    intersectZ = repmat(sensor.position(3), [size(newRays.origin, 1) 1]);
    intersectT = (intersectZ - newRays.origin(:, 3))./newRays.direction(:, 3);
    
    intersectPosition = newRays.origin + newRays.direction .* repmat(intersectT, [1 3]);
    imagePixel = [intersectPosition(:,2) intersectPosition(:, 1)]; 
    imagePixel = round(imagePixel* sensor.resolution(1)/sensor.size(1)    + repmat( sensor.resolution./2, [size(imagePixel,1) 1]));   %
    
    %make sure imagePixel is in range;
    imagePixel(imagePixel < 1) = 1;
    imagePixel = min(imagePixel, repmat(sensor.resolution, [size(imagePixel,1) 1]));
    
    %add a value to the intersection position
    for i = 1:size(croppedApertureSample.Y(:))
        sensor.image(imagePixel(i,1), imagePixel(i,2)) =  sensor.image(imagePixel(i,1), imagePixel(i,2)) + 1;  %sensor.image(imagePixel(:,1), imagePixel(:,2)) + 1;
    end
end

figure; imshow(sensor.image/ max(sensor.image(:)));
