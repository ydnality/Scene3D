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
clx
ieInit

%% Sum the wavefronts from a series of points in a slit aperture

% Assume single pinhole is at origin
% numSamples = 100000;
lambda    = 550;        % nm
binSize   = 2500;       % 25; % nm
numPixels = 1000;       % In the sensor

%apertureSample = 100;  % nm
%nApertureSamples = 1e3;

% We do the calculation one aperture sample at a time
% imagePlaneLocation = [lambda * 20 0];    % This many wavelengths from the aperture
imagePlaneDist     = lambda*2000;  %20
% imagePlane         = zeros(1,numPixels); % A 1D sensor
% ph = zeros(1, numPixels);               % Phase of ray at that sample
% distance = zeros(1, numPixels);         % From aperture sample to sensor sample

% intersectionPointR = zeros(1, numPixels);

% imagePlanePhase = zeros(1,numPixels);
% intersectionPoint = zeros(1, numPixels);

%we are using known incremental directions
%since all the rays that make it through are assumed to go through pinhole
%at the same phase, we can just start tracing from the pinhole

%aperturePoint = (rand(1,numApertureSamples) - .5) * 2 * apertureSize/2;

% Pick a set of sample points on the aperture
numApertureSamples = 100;   % Number of points in the slit aperture
apertureSize       = 2500;  % nm   25

aperturePoint = linspace(-.5,.5, numApertureSamples + 1) * apertureSize;

% These are the locations on the sensor
% A 2xN matrix (x,y) for the two rows
endLocations = linspace(-numPixels/2 * binSize, numPixels/2 * binSize, numPixels);
endLocations = [ones(1,numPixels)*imagePlaneDist;endLocations(:)'];

% For each aperture sample, calculate the distance and phase to each pixel
intensity = zeros(numPixels,1);

% For each aperture calculate the distance to all the pixels, keeping the
% phase information
for rPhase = 0 %:pi/256:(255 * pi/256)  % for rPhase = linspace(0, 2* pi, 14)
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

%% Sum the wavefronts from a series of points in a slit aperture 
% but this time, do a truly uniform spherical sample for huygens

% Assume single pinhole is at origin

% numSamples = 100000;
lambda    = 550;                % nm
% binSize   = 2500;   %25;                 % nm
% numPixels = 1000;                % In the sensor

binSize   = 2500;   %25;        % nm
numPixels = 250;                % In the sensor

% apertureSample = 100;           % nm
%nApertureSamples = 1e3;

% We do the calculation one aperture sample at a time
% imagePlaneLocation = [lambda * 20 0];    % This many wavelengths from the aperture
imagePlaneDist     = lambda*20000;  %20
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
% numApertureSamples = 100;   % Number of points in the slit aperture
% apertureSize = 2500;         % nm   25

numApertureSamples = 1000;   % Number of points in the slit aperture
apertureSize = 250000;         % nm   25

aperturePoint = linspace(-.5,.5, numApertureSamples + 1) * apertureSize;

% These are the locations on the sensor
% A 2xN matrix (x,y) for the two rows
endLocations = linspace(-numPixels/2 * binSize, numPixels/2 * binSize, numPixels);
endLocations = [ones(1,numPixels)*imagePlaneDist;endLocations(:)'];

% For rPhase = linspace(0, 2* pi, 14)
% For each aperture sample, calculate the distance and phase to each pixel
intensity = zeros(numPixels,1);

% For each aperture calculate the distance to all the pixels, keeping the
% phase information
for rPhase = 0 %:pi/256:(255 * pi/256)
    for apertureSample = 1:numApertureSamples
        
        % These are the locations on the aperture
        apLocations = [zeros(1,numPixels); ones(1, numPixels) * aperturePoint(apertureSample)];
        
        d = sqrt(diag((endLocations - apLocations)'*(endLocations - apLocations)));
        
        %sampling correction term
        dir = (endLocations - apLocations);
        normDir = normvec(dir, 'p', 2, 'dim', 1);
        
        %this factor accounts for the fact that we AREN'T sampling uniformly in theta.
        %theta = uniform
        %y = sin(theta)
        %use transformation of random variables to compute that fY(y) ~
        %1/(1-y^2). We take the inverse of this that is (1-y^2) since we
        %are sampling uniformly in y instead of theta.
        normalizeAmount = (sqrt(1 - normDir(2,:).^2))';
        
        thetaCapture = atan2(dir(2,:)+ binSize/2, dir(1,:) ) - atan2(dir(2,:)- binSize/2, dir(1,:) );
        
        % vcNewGraphWin; plot(d)
        % intensity = exp(2 * pi * 1i * (d/lambda)).*normalizeAmount + intensity;
        % intensity = exp(2 * pi * 1i * (d/lambda)).*thetaCapture'  + intensity;
        
        intensity = exp(2 * pi * 1i * (d/lambda))./d.^2  + intensity;
        %intensity = exp(2 * pi * 1i * (d/lambda)) + intensity;
    end
end

% Get the real intensity from the complex sum above
vcNewGraphWin;
intensity = abs(intensity).^2;

% Plot the intensity across the bin sizes.
sensorAxis = linspace(0, binSize * numPixels, numPixels);
plot(sensorAxis, intensity);

% Key observation: we don't always get an airy disk patern.  It depends
% heavily on the distance from the aperture, the aperture size, and the
% wavelength.
%
% See
% http://en.wikipedia.org/wiki/Fresnel_diffraction
% http://en.wikipedia.org/wiki/Fresnel_diffraction
% for more details.  We are seeing a transition between Fresnel (near
% field, and Fraunhofer

%% plot theoretical solution
theta = linspace(-pi/2, pi/2, 200);
a = 2500;
lambda = 550;
I = (sin(pi .* a .* sin(theta./lambda))./(pi * a .* sin(theta./lambda))).^2;

vcNewGraphWin;
plot (theta, I);

%% compare the 2 curves

thetaSens = atan((sensorAxis - sensorAxis(numPixels/2))/(imagePlaneDist));

vcNewGraphWin;
plot(thetaSens, intensity./max(intensity(:)));
hold on; plot(theta, I, 'r');
title('Analytical(red) vs. Huygens(blue) raytrace');

%why the error? the analytical solution uses a paraxial approximation...
%sin(theta) ~ tan(theta) for small angles

%% try in 3D instead

% assume single pinhole is at origin

% numSamples = 100000;
wave           = 400:150:700; %550;                % nm
binSize        = [7500 7500];   %25;                 % nm
numPixels      = [300 300];                % In the sensor
imagePlaneDist = 2000 * 550;  %20

% Pick a set of sample points on the aperture
numApertureSamples = [50 50];   % Number of points in the slit aperture
apertureSize = 2500;         % nm   25
numPixelsTot = numPixels(1) * numPixels(2);

%create aperture grid
ap1DX = linspace(-.5,.5, numApertureSamples(1) + 1) * apertureSize;
ap1DY = linspace(-.5,.5, numApertureSamples(2) + 1) * apertureSize;
[apXGrid, apYGrid] = meshgrid(ap1DX, ap1DY);
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
[endLGridX, endLGridY] = meshgrid(endLocations1DX, endLocations1DY);
endLGridXFlat = endLGridX(:);    %flatten the meshgrid
endLGridYFlat = endLGridY(:);
endLGridZFlat = ones(size(endLGridYFlat)) * imagePlaneDist;

intensity = zeros(numPixels(1), numPixels(2), length(wave));

oi = oiCreate;
oi = initDefaultSpectrum(oi);
oi = oiSet(oi, 'wave', wave);

countWave = 1;
%loop through all wavelenghts
for lambda = wave
    disp(lambda);
    
    intensityFlat = zeros(numPixels);
    intensityFlat = intensityFlat(:);
    %process 1 ray at the aperture at once
    for apertureSample = 1:numApertureSamplesTot
        %disp(apertureSample)
        % These are the locations on the aperture (we will take 1 at a time),
        % and repmat that.
        apLocationsX = ones(numPixelsTot, 1) * apXGridFlat(apertureSample);
        apLocationsY = ones(numPixelsTot, 1) * apYGridFlat(apertureSample);
        
        xDiff = endLGridXFlat - apLocationsX;
        yDiff = endLGridYFlat - apLocationsY;
        zDiff = endLGridZFlat; %aperture is assumed to be at Z = 0;
        
        %compute the distance from the apeture location
        d = sqrt(xDiff.^2 + yDiff.^2 + zDiff.^2);  %if we plug in a different d, then we can trace an entire lens
        
        %add the complex exponential contribution to the entire sensor and
        %sum
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

%% try using PSF camera

% Make a volume of point.  A little boring in this case
point = psCreate(0,0,-200);

% Read a lens file and create a lens
lensFileName = fullfile(s3dRootPath,'data', 'lens', 'dgauss.50mm.dat');

nSamples = 151; %501; %151;
apertureMiddleD = .5;   % mm
lens = lensC('apertureSample', [nSamples nSamples], ...
    'fileName', lensFileName, ...
    'apertureMiddleD', apertureMiddleD);

% Create a film (sensor)

% position - relative to center of final lens surface
% size - 'mm'
% wavelength samples
wave = lens.get('wave');

film = filmC('position', [0 0 1 ], ...
    'size', [1 1], ...
    'wave', wave);

% Ray trace the point to the film
camera = psfCameraC('lens',lens,'film',film,'point source',point{1});

% Sequence of events for estimating the PSF,
nLines = 100;
jitter = false;
%camera.estimatePSF(nLines,jitter);

%limits the entrance aperture so this can run faster
subsection = [-.025, -.025, .025, .025];
camera.computeRayDist(nLines,jitter,subsection);

% camera.estimatePSF(nLines,jitter,subsection);
% vcAddObject(camera.oiCreate); oiWindow;
%%  Huygens ray trace part - Calculate the sampling grid on the aperture
% This will become a function and better integrated

lambda         = 550;             %for now
binSize        = [40000 40000];   %25;                 % nm
numPixels      = [50 50];         % In the sensor
numPixelsTot   = numPixels(1) * numPixels(2);
imagePlaneDist = 16*1e6;         %100 * 10^6; %100mm
lensMode = false;

% Estimated that the width of the 1st zero of airy disk will be .0336
apXGridFlat = camera.rays.origin(:,1) * 10^6; %convert to nm from mm
apYGridFlat = camera.rays.origin(:,2) * 10^6;
apXGridFlat = apXGridFlat(~isnan(apXGridFlat));
apYGridFlat = apYGridFlat(~isnan(apYGridFlat));
numApertureSamplesTot = length(apXGridFlat);

vcNewGraphWin;
plot(apXGridFlat, apYGridFlat, 'o');  %plot the aperture samples
axis image

%% These are the locations on the sensor
endLocations1DX = linspace(-numPixels(1)/2 * binSize(1), numPixels(1)/2 * binSize(1), numPixels(1));
endLocations1DY = linspace(-numPixels(2)/2 * binSize(2), numPixels(2)/2 * binSize(2), numPixels(2));
[endLGridX, endLGridY] = meshgrid(endLocations1DX, endLocations1DY);
endLGridXFlat = endLGridX(:);    %flatten the meshgrid
endLGridYFlat = endLGridY(:);
endLGridZFlat = ones(size(endLGridYFlat)) * imagePlaneDist;

intensity = zeros(numPixels(1), numPixels(2), length(wave));
intensityFlat = zeros(numPixels);
intensityFlat = intensityFlat(:);

if (lensMode), initialD = camera.rays.distance(~isnan(camera.rays.distance));
else           initialD = zeros(numApertureSamplesTot, 1);
end

for apertureSample = 1:numApertureSamplesTot
    %     if(mod(apertureSample, 10) == 0)
    %         disp(apertureSample)
    %     end
    
    % These are the locations on the aperture (we will take 1 at a time),
    % and repmat that.
    apLocationsX = ones(numPixelsTot, 1) * apXGridFlat(apertureSample);
    apLocationsY = ones(numPixelsTot, 1) * apYGridFlat(apertureSample);
    
    xDiff = endLGridXFlat - apLocationsX;
    yDiff = endLGridYFlat - apLocationsY;
    zDiff = endLGridZFlat; %aperture is assumed to be at Z = 0;
    
    %compute the distance from the apeture location
    d = sqrt(xDiff.^2 + yDiff.^2 + zDiff.^2) + initialD(apertureSample) * 10^6;  %if we plug in a different d, then we can trace an entire lens
    
    %add the complex exponential contribution to the entire sensor and
    %sum
    intensityFlat = exp(2 * pi * 1i .* (d/lambda)) + intensityFlat;
end


%visualize - take note that we have a square root here to better visualize
intensityFlat = 1/lambda .* abs(intensityFlat).^2;
intensity1Wave = reshape(intensityFlat, size(intensity, 1), size(intensity, 2));

vcNewGraphWin;
imagesc(sqrt(intensity1Wave));
colormap('gray'); axis image

if(lensMode), lensModeText = 'Lens';
else          lensModeText = 'No Lens';
end
title([lensModeText ', '        int2str(imagePlaneDist/1e6) 'mm,'...
    ' apDiameter: '             num2str(apertureMiddleD) ...
    'mm, pointSourceLocation: ' num2str(-point{1}(3)) 'mm']);
xlabel([num2str(binSize(1) *numPixels(1)/1e6) 'mm']);


%% do it faster
tic
% Eventually want to limit this to several different passes because we will
% run out of memory...
xDiffMat = bsxfun(@minus, endLGridXFlat, apXGridFlat');
yDiffMat = bsxfun(@minus, endLGridYFlat, apYGridFlat');
zDiffMat = repmat(endLGridZFlat, [1, numApertureSamplesTot]);
expD = exp(2 * pi * 1i .* (sqrt(xDiffMat.^2 + yDiffMat.^2 + zDiffMat.^2)/lambda));
intensityFlat = sum(expD, 2) + intensityFlat;

intensityFlat = 1/lambda .* abs(intensityFlat).^2;
intensity1Wave = reshape(intensityFlat, size(intensity, 1), size(intensity, 2));

vcNewGraphWin;
imagesc(sqrt(intensity1Wave));
colormap('gray'); axis image

if(lensMode), lensModeText = 'Lens';
else          lensModeText = 'No Lens';
end
title([lensModeText ', '        int2str(imagePlaneDist/1e6) 'mm,'...
    ' apDiameter: '             num2str(apertureMiddleD) ...
    'mm, pointSourceLocation: ' num2str(-point{1}(3)) 'mm']);
xlabel([num2str(binSize(1) *numPixels(1)/1e6) 'mm']);

toc
%% do it faster but better memory management

tic

intensityFlat = zeros(numPixels);
intensityFlat = intensityFlat(:);
jobInterval   = 5000;
numJobs       = ceil(numApertureSamplesTot/jobInterval);

%split aperture into segments and do bsxfun on that, then combine later.
%This can be parallelized later
for job = 1:numJobs
    
    apXGridFlatCur = apXGridFlat((job-1) * jobInterval + 1:min(job * jobInterval, numApertureSamplesTot));
    apYGridFlatCur = apYGridFlat((job-1) * jobInterval + 1:min(job * jobInterval, numApertureSamplesTot));
    
    xDiffMat = bsxfun(@minus, endLGridXFlatCur, apXGridFlatCur');
    yDiffMat = bsxfun(@minus, endLGridYFlatCur, apYGridFlatCur');
    zDiffMat = repmat(endLGridZFlat, [1, numApertureSamplesTot]);
    
    expD = exp(2 * pi * 1i .* (sqrt(xDiffMat.^2 + yDiffMat.^2 + zDiffMat.^2)/lambda));
    intensityFlat = sum(expD, 2) + intensityFlat;
end

intensityFlat = abs(intensityFlat);

toc

%% test the new infrastructure with PSFCamera and estimatePSF

% Make a volume of point.  A little boring in this case
point = psCreate(0,0,-200);

% Read a lens file and create a lens
lensFileName = fullfile(s3dRootPath,'data', 'lens', 'dgauss.50mm.dat');

nSamples = 501; %501; %151;
apertureMiddleD = .1;  %.5;   % mm
lens = lensC('apertureSample', [nSamples nSamples], ...
    'fileName', lensFileName, ...
    'apertureMiddleD', apertureMiddleD, ...
    'diffractionEnabled', true);

% Create a film (sensor)
%
% position - relative to center of final lens surface
% size - 'mm'
% wavelength samples
lens.set('wave', 400:150:700);
wave = lens.get('wave');

%put it 16 mm away
film = filmC('position', [0 0 16 ], ...
    'resolution', [50 50 1], ...
    'size', [2/2 2/2], ...
    'wave', wave);

% Ray trace the point to the film
camera = psfCameraC('lens',lens,'film',film,'point source',point{1});

% Sequence of events for estimating the PSF,
nLines = 100;
jitter = false;

%limits the entrance aperture so this can run faster
%subsection = [-.03, -.03, .03, .03];  % for aperture .5
subsection = [-.03/5, -.03/5, .03/5, .03/5];
method = 'huygens';
rtType = 'realistic';
camera.estimatePSF(nLines, jitter, subsection, method, rtType);
oi = camera.oiCreate(); 
vcAddObject(oi); oiWindow;


%% Produce diffraction pattern for the realistic lens test case

% Make a volume of point.  A little boring in this case
point = psCreate(0,0,-9e7);

% Read a lens file and create a lens
% lensFileName = fullfile(s3dRootPath,'data', 'lens', 'dgauss.50mm.dat');
lensFileName = fullfile(s3dRootPath,'data', 'lens', '2ElLens.dat');

% nSamples = 4001; %501; %151;
% apertureMiddleD = 2.2727;  %.5;   % mm

nSamples        = 101;     % 501; %151;
apertureMiddleD = .22727;  %.5;   % mm
lens = lensC('apertureSample', [nSamples nSamples], ...
    'fileName', lensFileName, ...
    'apertureMiddleD', apertureMiddleD, ...
    'diffractionEnabled', true);

% Create a film (sensor)
% position - relative to center of final lens surface
% size - 'mm'
% wavelength samples
lens.set('wave', 400:100:700);
wave = lens.get('wave');

%put it 16 mm away
film = filmC('position', [0 0 48.5 ], ...
    'resolution', [25 25 1], ...
    'size', [2/sqrt(2) 2/sqrt(2)], ...
    'wave', wave);

% Ray trace the point to the film
camera = psfCameraC('lens',lens,'film',film,'point source',point{1});

% Sequence of events for estimating the PSF,
nLines = 100;
jitter = false;
%camera.estimatePSF(nLines,jitter);

%limits the entrance aperture so this can run faster
%subsection = [-.03, -.03, .03, .03];  % for aperture .5
subsection = [-.25, -.25, .25, .25];

method = 'huygens';
rtType = 'realistic';
camera.estimatePSF(nLines,jitter,subsection, method, rtType);
oi = camera.oiCreate(); 
vcAddObject(oi); oiWindow;

%% Produce diffraction pattern using an ideal lens

% Make a volume of point.  A little boring in this case
point = psCreate(0,0,-1e17);

% Read a lens file and create a lens
%lensFileName = fullfile(s3dRootPath,'data', 'lens', 'dgauss.50mm.dat');
lensFileName = fullfile(s3dRootPath,'data', 'lens', '2ElLens.dat');

% Not sure what AL means by work, not work, here.
% nSamples = 5001; %501; %151;
% apertureMiddleD = 2.2727;  %.5;   % mm   %DID NOT WORK - WHY?

% nSamples = 1001; %501; %151;
% apertureMiddleD = .44;  %.5;   % mm     %did not work - WHY    1/sqrt(2)
%size sensor

%nSamples = 251; %501; %151;
%apertureMiddleD = .11;  %.5;   % mm    %WORKS BRILLIANTLY

nSamples        = 501;     %501; %151;
apertureMiddleD = .22727;  %.5;   % mm    %WORKS BRILLIANTLY
lens = lensC('apertureSample', [nSamples nSamples], ...
    'fileName', lensFileName, ...
    'apertureMiddleD', apertureMiddleD, ...
    'diffractionEnabled', true);

% position - relative to center of final lens surface
% size - 'mm'
% wavelength samples
lens.set('wave', 400:100:700);
wave = lens.get('wave');

% Create a film (sensor)
% put it 16 mm away
film = filmC('position', [0 0 50], ...
    'resolution', [100 100 1], ...
    'size', [1/sqrt(2) 1/sqrt(2)], ...
    'wave', wave);

% Ray trace the point to the film
camera = psfCameraC('lens',lens,'film',film,'point source',point{1});

% Sequence of events for estimating the PSF,
nLines = 100;
jitter = false;

%Setting subsection limits the entrance aperture so this can run faster
%
%subsection = [-.03, -.03, .03, .03];  % for aperture .5
%subsection = [-.25, -.25, .25, .25];
%subsection = [-1 -1, 1, 1];
% When empty, the whole lens is used.
subsection = [];

% Diffraction method
% method = 'HURB';
method = 'huygens';

% Lens type for ray trace is ideal
rtType = 'ideal';
camera.estimatePSF(nLines,jitter,subsection, method, rtType);
oi = camera.oiCreate();
vcAddObject(oi); oiWindow;
%% END

%% observations

%aperture szie: .1mm
%sensor size:.25  mm
%sensor position (10*10^6 nm)
%barely any difference.  the distance tracking one seems more "Focused"
%which could make sense.
%pattern does not look like airy disk.  This might be ok since we are in
%near field diffraction


