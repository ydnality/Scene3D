%% test simple ray-wave-tracing 
%
%  We trace waves from a single slit aperture to a point on the sensor.  As
%  each ray emerges from the slit aperture towards the sensor, it gives
%  rise to a spherical wavefront.  We sample rays in various directions on
%  the sphere and keep track of the distance they travel from the aperture
%  to the sensor.
%
%  We know the wavelength of the light, and we assume the source is
%  coherent (all in the same phase when arriving at the aperture.  So the
%  amplitude at the sensor is given by a cosine function.  This gives rise
%  to an interference pattern at the sensor.
%
%
%
%
% AL
%%
ieInit

%% Sum the wavefronts from a series of points in a slit aperture

%assume single pinhole is at origin

% numSamples = 100000;
lambda    = 550;                % nm
binSize   = 2500;   %25;                 % nm
numPixels = 1000;                % In the sensor
%apertureSample = 100;           % nm
%nApertureSamples = 1e3;

% We do the calculation one aperture sample at a time
% imagePlaneLocation = [lambda * 20 0];    % This many wavelengths from the aperture
imagePlaneDist     = lambda*2000;  %20
% imagePlane         = zeros(1,numPixels); % A 1D sensor
ph = zeros(1, numPixels);               % Phase of ray at that sample
distance = zeros(1, numPixels);         % From aperture sample to sensor sample

% intersectionPointR = zeros(1, numPixels);

% imagePlanePhase = zeros(1,numPixels);   
% intersectionPoint = zeros(1, numPixels);

%we are using known incremental directions
%since all the rays that make it through are assumed to go through pinhole
%at the same phase, we can just start tracing from the pinhole

%aperturePoint = (rand(1,numApertureSamples) - .5) * 2 * apertureSize/2;

% Pick a set of sample points on the aperture
numApertureSamples = 100;   % Number of points in the slit aperture
apertureSize = 2500;         % nm   25
aperturePoint = linspace(-.5,.5, numApertureSamples + 1) * apertureSize;

% These are the locations on the sensor
% A 2xN matrix (x,y) for the two rows
endLocations = linspace(-numPixels/2 * binSize, numPixels/2 * binSize, numPixels);
endLocations = [ones(1,numPixels)*imagePlaneDist;endLocations(:)'];

%for rPhase = linspace(0, 2* pi, 14)
% For each aperture sample, calculate the distance and phase to each pixel
intensity = zeros(numPixels,1);

% For each aperture calculate the distance to all the pixels, keeping the
% phase information
for rPhase = 0 %:pi/256:(255 * pi/256)
    for apertureSample = 1:numApertureSamples

        % These are the locations on the aperture
        apLocations = [zeros(1,numPixels); ones(1, numPixels) * aperturePoint(apertureSample)];

        d = sqrt(diag((endLocations - apLocations)'*(endLocations - apLocations)));
        % vcNewGraphWin; plot(d)
        intensity = exp(2 * pi * 1i * (d/lambda)) + intensity;
    end
end
% Get the real intensity
vcNewGraphWin;
intensity = abs(intensity).^2;

sensorAxis = linspace(0, binSize * numPixels, numPixels);
plot(sensorAxis, intensity);

%% plot theoretical solution
theta = linspace(-pi/2, pi/2, 200);
a = 2500;
lambda = 550;
I = (sin(pi .* a .* sin(theta./lambda))./(pi * a .* sin(theta./lambda))).^2;
figure; 
plot (theta, I);

%% compare the 2 curves
thetaSens = atan((sensorAxis - sensorAxis(numPixels/2))/(imagePlaneDist));
vcNewGraphWin;
plot(thetaSens, intensity./max(intensity(:)));
hold on;
plot(theta, I, 'r');

title('Analytical(red) vs. Huygens(blue) raytrace');

%why the error? the analytical solution uses a paraxial approximation...
%sin(theta) ~ tan(theta) for small angles

%% try in 3D instead

%assume single pinhole is at origin

% numSamples = 100000;
wave    = 400:100:700;%550;                % nm
binSize   = [7500 7500];   %25;                 % nm
numPixels = [300 300];                % In the sensor
imagePlaneDist     = 2000 * 550;  %20

% Pick a set of sample points on the aperture
numApertureSamples = [50 50];   % Number of points in the slit aperture
apertureSize = 2500;         % nm   25
numPixelsTot = numPixels(1) * numPixels(2);

%create aperture grid
ap1DX = linspace(-.5,.5, numApertureSamples(1) + 1) * apertureSize;
ap1DY = linspace(-.5,.5, numApertureSamples(2) + 1) * apertureSize;
[apXGrid apYGrid] = meshgrid(ap1DX, ap1DY);
apXGridFlat = apXGrid(:);   %flatten to 1D for simplicity
apYGridFlat = apYGrid(:);
%reduce aperture grid by applying a circular aperture.  Try different
%apertures!!
withinAperture = (apXGridFlat.^2 + apYGridFlat.^2) < (apertureSize/2)^2;
%withinAperture = logical(ones(size(apXGridFlat)));
apXGridFlat = apXGridFlat(withinAperture);
apYGridFlat = apYGridFlat(withinAperture);
numApertureSamplesTot = (sum(withinAperture));

%create sensor grid
% These are the locations on the sensor
endLocations1DX = linspace(-numPixels(1)/2 * binSize(1), numPixels(1)/2 * binSize(1), numPixels(1));
endLocations1DY = linspace(-numPixels(2)/2 * binSize(2), numPixels(2)/2 * binSize(2), numPixels(2));
[endLGridX endLGridY] = meshgrid(endLocations1DX, endLocations1DY);
endLGridXFlat = endLGridX(:);    %flatten the meshgrid
endLGridYFlat = endLGridY(:);
endLGridZFlat = ones(size(endLGridYFlat)) * imagePlaneDist;

intensity = zeros(numPixels(1), numPixels(2), length(wave));

oi = oiCreate;
oi = initDefaultSpectrum(oi);
oi = oiSet(oi, 'wave', wave);

countWave = 1;
for lambda = wave
    disp(lambda);
    
    intensityFlat = zeros(numPixels);
    intensityFlat = intensityFlat(:);
    for apertureSample = 1:numApertureSamplesTot
        %disp(apertureSample)
        % These are the locations on the aperture (we will take 1 at a time),
        % and repmat that.
        apLocationsX = ones(numPixelsTot, 1) * apXGridFlat(apertureSample);
        apLocationsY = ones(numPixelsTot, 1) * apYGridFlat(apertureSample);

        xDiff = endLGridXFlat - apLocationsX;
        yDiff = endLGridYFlat - apLocationsY;
        zDiff = endLGridZFlat; %aperture is assumed to be at Z = 0;

        d = sqrt(xDiff.^2 + yDiff.^2 + zDiff.^2);

        intensityFlat = exp(2 * pi * 1i .* (d/lambda)) + intensityFlat;
    end
    
    intensityFlat = 1/lambda .* abs(intensityFlat).^2;
    intensity1Wave = reshape(intensityFlat, size(intensity, 1), size(intensity, 2));
    
    % sensorAxisX = linspace(0, binSize(1) * numPixels(1), numPixels(1));
    % sensorAxisY = linspace(0, binSize(2) * numPixels(2), numPixels(2));
    % 
    % vcNewGraphWin; imagesc(intensity);
    % vcNewGraphWin;  surf(sensorAxisX(1:2:300), sensorAxisY(1:2:300), intensity(1:2:300, 1:2:300), 'LineStyle', 'none');
    %[sensorAxisXGrid sensorAxisYGrid] = meshgrid(sensorAxisX, sensorAxisY);
    %plot(sensorAxis, intensity);
    
    intensity(:,:, countWave) = intensity1Wave;
    countWave = countWave + 1;
end

oi = oiSet(oi, 'photons', intensity);
vcAddObject(oi); oiWindow;