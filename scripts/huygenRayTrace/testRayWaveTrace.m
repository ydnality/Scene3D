%% test simple ray-wave-tracing 

%% The first example is a simple point source, through a single pinhole slit.  We are tracing in the forward direction in this case!  Only 2D for now

%assume point source is at (-10, 0)
%assume single pinhole is at origin
%assume imaging plane is at (10, 0)

initialIntensity = 10;
numSamples = 300000;
numApertureSamples = 80;
lambda = 550;  %nm
binSize = 1/5; %nm
apertureSize = 10;
numPixels = 5000;

imagePlaneLocation = [lambda * 20 0];
imagePlane = zeros(1,numPixels);
phase = zeros(1, numSamples);
distance = zeros(1, numSamples);
intersectionPointR = zeros(1, numSamples);
intersectionPoint = zeros(1, numSamples);
intensity = zeros(1, numSamples);

%we are using known incremental directions
%since all the rays that make it through are assumed to go through pinhole
%at the same phase, we can just start tracing from the pinhole
for apertureSample = 1:numApertureSamples
    
    randSamples = rand(1,numSamples);
    for i = 1:numSamples
        theta = (randSamples(i) - .5) * pi/2;
        rayVector = [cos(theta) sin(theta)];

        %compute interesction point with plane, and put intensity onto the
        %plane
        distance(i) = imagePlaneLocation(1)/rayVector(1);
        intersectionPoint(i) = rayVector(2)*distance(i) + (apertureSample - numApertureSamples/2) * apertureSize/2 ;

        %condition the intersectionPoint so it fits in the 1D image plane;
        phase(i) = (mod(distance(i), lambda)/lambda);
        intersectionPointR(i) = round(intersectionPoint(i) * binSize + size(imagePlane,2)/2);
%         if (intersectionPointR(i) > size(imagePlane,2))
%             intersectionPointR(i) = size(imagePlane,2);
%         elseif(intersectionPointR(i) < 1)
%             intersectionPointR(i) = 1;   
%         end

        %if ray is still within bounds, record it
        if(intersectionPointR(i) <= size(imagePlane,2) && intersectionPointR(i) >=1)
            intensity(i) = cos(2*pi * phase(i));
            %imagePlane(intersectionPointR(i)) = imagePlane(intersectionPointR(i)) + intensity(i);
            imagePlane(intersectionPointR(i)) = imagePlane(intersectionPointR(i)) + intensity(i);
        end
    end
end

sensorAxis = linspace(0, binSize * numPixels, numPixels);
figure; plot(sensorAxis, abs(imagePlane.^2));

%% try in 3D instead

