%% Ray-tracing for an ideal lens
%
% What does this do?
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

% Declare point sources in world space.  The camera is usually at [0 0 0],
% and pointing towards -z.  We are using a right-handed coordinate system.
% pointSources = [0 0 -100;
%                 0 50 -100;
%                 0 -50 -100;
%                 50 0 -100;
%                 50 50 -100;
%                 50 -50 -100;
%                 -50 0 -100;
%                 -50 50 -100;
%                 -50 -50 -100;];

%%  Create the array of point sources.
%

[XGrid YGrid] = meshgrid(-4000:1000:4000,-4000:1000:4000);

% Points are on a plane in 3D
% Not sure about the units.
pointSources = [XGrid(:) YGrid(:) ones(size(XGrid(:))) * -20000];

%% Sensor properties

% This should be integrated with the ISET sensor idea.  Or given a
% different name.  Or both.
sensor.size = [48 48]; %in mm
sensor.position = [0 0 50]; 
sensor.resolution = [200 200];
sensor.image = zeros(sensor.resolution);
sensor.focalLength = 50; %in mm

apertureRadius =5; % in mm

% this position should NOT change
lensCenterPosition = [0 0 0];

%loop through all point sources
for curInd = 1:size(pointSources, 1);
    
    % This calculation happens a lot ... we should functionalize it.
    curPointSource = pointSources(curInd, :);
    
    %trace ray from point source to lens center, to image.  This helps
    %determine the point of focus
    centerRay.origin = curPointSource;
    centerRay.direction = lensCenterPosition - centerRay.origin;
    centerRay.direction = centerRay.direction./norm(centerRay.direction);
    
    %calculate the in-focus plane using thin lens equation
    inFocusDistance = 1/(1/sensor.focalLength - -1/curPointSource(3));
    
    % What is this?
    inFocusT = (inFocusDistance - centerRay.origin(3))/centerRay.direction(3);
    inFocusPosition = centerRay.origin + inFocusT .* centerRay.direction;
    
    %Create a set of circular aperture positions and uniformly sample the aperture
    %everything.  Calculations will be done in vector form for speed
    % The code here is a good candidate for a function (BW).
    [apertureSample.X, apertureSample.Y] = meshgrid(linspace(-1, 1, 90),linspace(-1, 1, 90)); %adjust this if needed
    
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
       
    %when intersecting ideal lens, change the direction to intersect the
    %inFocusPosition, and update the origin
    lensIntersectT = (lensCenterPosition(3) - rays.origin(:,3))./ rays.direction(:,3);
    lensIntersectPosition = rays.origin +  repmat(lensIntersectT, [1 3]) .* rays.direction;
    %calculate new direction
    newRays.origin = lensIntersectPosition;
    newRays.direction = repmat(inFocusPosition , [size(newRays.origin,1) 1 ]) - newRays.origin;
    
    %calculate intersection point at sensor
    intersectZ = repmat(sensor.position(3), [size(newRays.origin, 1) 1]);
    intersectT = (intersectZ - newRays.origin(:, 3))./newRays.direction(:, 3);
    
    intersectPosition = newRays.origin + newRays.direction .* repmat(intersectT, [1 3]);
    imagePixel = [intersectPosition(:,2) intersectPosition(:, 1)]; 
    imagePixel = round(imagePixel * sensor.resolution(1)/sensor.size(1)  + repmat( sensor.resolution./2, [size(imagePixel,1) 1]));
    
    %make sure imagePixel is in range;
    imagePixel(imagePixel < 1) = 1;
    imagePixel = min(imagePixel, repmat(sensor.resolution, [size(imagePixel,1) 1]));
    
    %add a value to the intersection position
    for i = 1:size(croppedApertureSample.Y(:))
        sensor.image(imagePixel(i,1), imagePixel(i,2)) =  sensor.image(imagePixel(i,1), imagePixel(i,2)) + 1;  %sensor.image(imagePixel(:,1), imagePixel(:,2)) + 1;
    end
end

%%
vcNewGraphWin; imshow(sensor.image/ max(sensor.image(:)));

%%