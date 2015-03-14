%% test simple ray-wave-tracing 

%% The first example is a simple point source, through a single pinhole slit.  We are tracing in the forward direction in this case!  Only 2D for now

%assume point source is at (-10, 0)
%assume single pinhole is at origin
%assume imaging plane is at (10, 0)

initialIntensity = 10;
numSamples = 100000;
numApertureSamples = 100;
lambda = 550;  %nm
binSize = 128; %nm   %1/2
apertureSize =4000; %nm
numPixels = 625;

imagePlaneLocation = [lambda * 200 0];   %lambda * 20
imagePlane = zeros(1,numPixels);
imagePlanePhase = zeros(1,numPixels);
phase = zeros(1, numSamples);
distance = zeros(1, numSamples);
intersectionPointR = zeros(1, numSamples);
intersectionPoint = zeros(1, numSamples);
intensity = zeros(1, numSamples);
%we are using known incremental directions
%since all the rays that make it through are assumed to go through pinhole
%at the same phase, we can just start tracing from the pinhole

aperturePoint = (rand(1,numApertureSamples) - .5) * 2 * apertureSize/2;
aperturePoint = linspace(-.5,.5, numApertureSamples) * 2 * apertureSize/2;

rPhase = 0;
    for apertureSample = 1:numApertureSamples

        randSamples = rand(1,numSamples);

        for i = 1:numSamples
            
            %naive theta - shoot them at a pi/2 cone
            theta = (randSamples(i) - .5) * pi/2;

            rayVector = [cos(theta) sin(theta)];

   
            %compute interesction point with plane, and put intensity onto the
            %plane
            distance(i) = imagePlaneLocation(1)/rayVector(1);
            %distance(i) = imagePlaneLocation(1);  %spherical sensor

            %intersectionPoint(i) = rayVector(2)*distance(i) + (apertureSample - numApertureSamples/2) * apertureSize/2 ;

            intersectionPoint(i) = rayVector(2)*distance(i) + aperturePoint(apertureSample);
            
            %condition the intersectionPoint so it fits in the 1D image plane;
            phase(i) = (mod(distance(i), lambda)/lambda);

            intersectionPointR(i) = round(intersectionPoint(i)/binSize + size(imagePlane,2)/2);

            %if ray is still within bounds, record it
            if(intersectionPointR(i) <= size(imagePlane,2) && intersectionPointR(i) >=1)
                intensity(i) = cos(2*pi * -phase(i) + rPhase);

                %imagePlane(intersectionPointR(i)) = imagePlane(intersectionPointR(i)) + intensity(i);

                imagePlane(intersectionPointR(i)) = imagePlane(intersectionPointR(i)) + intensity(i);
            end
            count = count + 1;
        end
    end


sensorAxis = linspace(0, binSize * numPixels, numPixels);
figure; plot(sensorAxis, abs(imagePlane.^2));


%% faster version of above

%assume point source is at (-10, 0)
%assume single pinhole is at origin
%assume imaging plane is at (10, 0)

numSamples = 100000;
numApertureSamples = 5000;
lambda = 550;  %nm
binSize = 20; %nm   %1/2
apertureSize =2000; %nm
numPixels = 625;

imagePlaneLocation = [lambda * 20 0];   %lambda * 20
imagePlane = zeros(1,numPixels);
imagePlanePhase = zeros(1,numPixels);
phase = zeros(1, numSamples);
distance = zeros(1, numSamples);
intersectionPointR = zeros(1, numSamples);
intersectionPoint = zeros(1, numSamples);
intensity = zeros(1, numSamples);
%we are using known incremental directions
%since all the rays that make it through are assumed to go through pinhole
%at the same phase, we can just start tracing from the pinhole

%aperturePoint = (rand(1,numApertureSamples) - .5) * 2 * apertureSize/2;
aperturePoint = linspace(-.5,.5, numApertureSamples) * 2 * apertureSize/2;


numSamples = numPixels;  %smarter way to trace rays... 

%for  rPhase = linspace(-pi, pi, 1)
rPhase = 0;
%for rPhase = linspace(0, 2* pi, 14)
    for apertureSample = 1:numApertureSamples

       % randSamples = rand(1,numSamples);
        
        %smarter theta - only shoot rays AT the sensor
        endLocations = linspace(-numPixels/2 * binSize, numPixels/2 * binSize, numPixels);
        
        %figure out thetas given this
        dir = [ones(1, numPixels) * imagePlaneLocation(1); endLocations] - [zeros(1,numPixels); ones(1, numPixels) * aperturePoint(apertureSample)];
        dir = dir ./ repmat(sqrt(sum(dir.*dir, 1)), [2 1]); %normalize
        
        count = 1;
        for i = 1:numSamples
            rayVector = dir(:, i);  %new stuff
            distance(i) = imagePlaneLocation(1)/rayVector(1);

            %phase(i) = (mod(distance(i), lambda)/lambda);
            phase(i) = distance(i)/lambda;
            intersectionPointR(i) = count;  %endLocations(i);
            
            %if ray is still within bounds, record it
           % if(intersectionPointR(i) <= size(imagePlane,2) && intersectionPointR(i) >=1)
                intensity(i) = cos(2*pi * phase(i) + rPhase);
                imagePlane(intersectionPointR(i)) = imagePlane(intersectionPointR(i)) + intensity(i);
            %end
            count = count + 1;
        end
    end
    %rPhase
%end

% theta = atan(endLocations./imagePlaneLocation(1)) * 180./pi; 
% figure; plot(theta, endLocations);

sensorAxis = linspace(0, binSize * numPixels, numPixels);
figure; plot(sensorAxis, abs(imagePlane.^2));

%% make it faster
theta = linspace(-pi/2, pi/2, 200);
a = 2000;
lambda = 550;
I = (sin(pi .* a .* sin(theta./lambda))./(pi * a .* sin(theta./lambda))).^2;
figure; 
plot (theta, I);
% plot theoretical solutiongoo


%% try in 3D instead

