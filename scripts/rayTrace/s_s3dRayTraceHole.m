%% ray-tracing for an aperture only


% declare point sources in world space.  The camera is usually at [0 0 0],
% and pointing towards -z.  We are using a right-handed coordinate system.
pointSources = [0 0 -10;
                0 5 -10;
                0 -5 -10;
                5 0 -10;
                5 5 -10;
                5 -5 -10;
                -5 0 -10;
                -5 5 -10;
                -5 -5 -10;];
 

%sensor properties
sensor.size = [24 24]; %in mm
sensor.position = [0 0 100]; 
sensor.resolution = [200 200];
sensor.image = zeros(sensor.resolution);
% sensor.focalLength = 50; %in mm

apertureRadius =1; % in mm

% this position should NOT change
lensCenterPosition = [0 0 0];

%loop through all point sources
for curInd = 1:size(pointSources, 1);
    
    curPointSource = pointSources(curInd, :);
    
    %trace ray from point source to lens center, to image.  This helps
    %determine the point of focus
    centerRay.origin = curPointSource;
    centerRay.direction = lensCenterPosition - centerRay.origin;
    centerRay.direction = centerRay.direction./norm(centerRay.direction);
    
    %loop through aperture positions and uniformly sample the aperture
    %everything is done in vector form for speed
    [apertureSample.X, apertureSample.Y] = meshgrid(linspace(-1, 1, 90),linspace(-1, 1, 90)); %adjust this if needed
    
    %assume a circular aperture, and make a mask that is 1 when the pixel
    %is within a circle of radius .5
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
    
    %calculate intersection point at sensor
    intersectZ = repmat(sensor.position(3), [size(rays.origin, 1) 1]);
    intersectT = (intersectZ - rays.origin(:, 3))./rays.direction(:, 3);
    
    intersectPosition = rays.origin + rays.direction .* repmat(intersectT, [1 3]);
    imagePixel = [intersectPosition(:,2) intersectPosition(:, 1)]; 
    imagePixel = round(imagePixel + repmat( sensor.resolution./2, [size(imagePixel,1) 1]));
    
    %make sure imagePixel is in range;
    imagePixel(imagePixel < 1) = 1;
    imagePixel = min(imagePixel, repmat(sensor.resolution, [size(imagePixel,1) 1]));
    
    %add a value to the intersection position
    for i = 1:size(croppedApertureSample.Y(:))
        sensor.image(imagePixel(i,1), imagePixel(i,2)) =  sensor.image(imagePixel(i,1), imagePixel(i,2)) + 1;  %sensor.image(imagePixel(:,1), imagePixel(:,2)) + 1;
    end
end

figure; imshow(sensor.image/ max(sensor.image(:)));
