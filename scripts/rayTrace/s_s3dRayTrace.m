%% Point source ray-tracing prototype
%
% What does this do?
%   Illustrate creating a ray structure - TODO, make a rayCreate/Set/Get?
%   Illustrate geometric calculation with rays
%
% (AL) Vistalab, 2013 

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
 
%% pinhole camera test

%sensor properties
sensor.size = [24 24]; %in mm
sensor.position = [0 0 100]; 
sensor.resolution = [200 200];
sensor.image = zeros(sensor.resolution);

% simple pinhole test
pinholePosition = [0 0 0];

%loop through all point sources
for curInd = 1:size(pointSources, 1);
    
    curPointSource = pointSources(curInd, :);
    
    %trace ray from point source to pinhole, to image.  only 1 sample is
    %needed in the pinhole case
    ray.origin = curPointSource;
    ray.direction = pinholePosition - ray.origin;
    ray.direction = ray.direction./norm(ray.direction);
    
    %calculate intersection point at sensor
    intersectZ = sensor.position(3);
    intersectT = (intersectZ - ray.origin(3))/ray.direction(3);
    
    intersectPosition = ray.origin + ray.direction .* intersectT;
    imagePixel = [intersectPosition(2) intersectPosition(1)]; 
    imagePixel = round(imagePixel + sensor.resolution./2);
    
    %make sure imagePixel is in range;
    imagePixel(imagePixel < 1) = 1;
    imagePixel = min(imagePixel, sensor.resolution);
    
    %add a value to the intersection position
    sensor.image(imagePixel(1), imagePixel(2)) = sensor.image(imagePixel(1), imagePixel(2)) + 1;
end

%%
vcNewGraphWin; imshow(sensor.image);


%% End


