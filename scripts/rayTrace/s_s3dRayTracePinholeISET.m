%% point source ray-tracing prototype
%
%  Try to move this stuff into ISET land so we get the other benefits, such
%  as spectral management of the sensor.  Think about this.
%
% BW/AL Vistalab 2014

%% Create point sources

% Should we be able to do this as a scene structure?

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
 
%% Pinhole camera test

% simple pinhole test
pinholePosition = [0 0 0];

% This means distance behind the pinhole.  In general it is the sensor
% position with respect to the focal length of the lens?
sensor.position = [0 0 100];     % This is position in 3 space, I think.

%sensor properties
sensor.size = [24 24]; %in mm

sensor.resolution = [200 200];   % What is this?  Number of pixels, I guess
sensor.image = zeros(sensor.resolution);

%loop through all point sources
for curInd = 1:size(pointSources, 1);
    
    % Here is the first point position
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
    
    %add a irradiance value to the intersection position on the sensor
    %surface. 
    sensor.image(imagePixel(1), imagePixel(2)) = sensor.image(imagePixel(1), imagePixel(2)) + 1;
end

%%
vcNewGraphWin;
imshow(sensor.image);

%% End

