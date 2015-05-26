%% 2014 OSA Conference
%
% This script uses the point test file, blur it in 2 different ways.  The
% first way is to use the new ray tracing method which uses Heisenburg
% Uncertainty Ray Bending (HURB).  The second way is the classical way,
% using theoretical PSF's.
%
%  NEEDS TO BE FIXED FOR FINDING THE DATA
%
% AL

%% Initialize ISET
ieInit;

%% Specify HURB ray tracing location and specification
% chdir(fullfile(s3dRootPath, 'papers', '2014-OSA'));
% sampleArray = cell(1, 1);
% 
% sampleArray{1}.rayTraceFile = 'PSFCenter_50mm_2m_f22_n401.mat'%'25mm_1m_65res.pbrt.mat' %'rayTrace25mm32res.mat' 
% sampleArray{1}.focalLength = 50
% sampleArray{1}.apertureDiameter = 2.2727
% sampleArray{1}.filmDistance = 51.2821	
% sampleArray{1}.targetDistance = 2
% 

%% Produce HURB results
% Make a point source (approximately infinity)
point = psCreate(0,0,-1000000000);

% Read a lens file and create a lens
%lensFileName = fullfile(s3dRootPath,'data', 'lens', 'dgauss.50mm.dat');
lensFileName = fullfile(s3dRootPath,'data', 'lens', '2ElLens.dat');

nSamples = 2001; 
apertureMiddleD = .11;  %.5;   % mm    %WORKS BRILLIANTLY

lens = lensC('apertureSample', [nSamples nSamples], ...
    'fileName', lensFileName, ...
    'apertureMiddleD', apertureMiddleD, ...
    'diffractionEnabled', true);

% Create a film (sensor) 
% position - relative to center of final lens surface
% size - 'mm'
% wavelength samples
lens.set('wave', [400:10:700 ]);

wave = lens.get('wave');

%put it 16 mm away
film = filmC('position', [0 0 50], ...
    'resolution', [300 300 1], ...
    'size', [2/sqrt(2) 2/sqrt(2)], ...
    'wave', wave);

% Create a camera out of lens, film ,and point source
camera = psfCameraC('lens',lens,'film',film,'point source',point{1});

% Sequence of events for estimating the PSF, 
nLines = 100;
jitter = false;
%camera.estimatePSF(nLines,jitter);

%limits the entrance aperture so this can run faster
subsection = [];

method = 'HURB';
%method = 'huygens';
rtType = 'ideal';
camera.estimatePSF(nLines,jitter,subsection, method, rtType);

oiHURB = camera.oiCreate(); vcAddObject(oiHURB); oiWindow;
    

%% Produce Huygens-Fresnel results (this section takes a long time to run)

% Make a point source (approximately infinity)
point = psCreate(0,0,-1000000000000000000);

% Read a lens file and create a lens
%lensFileName = fullfile(s3dRootPath,'data', 'lens', 'dgauss.50mm.dat');
lensFileName = fullfile(s3dRootPath,'data', 'lens', '2ElLens.dat');

nSamples = 401; %501; %151;
apertureMiddleD = .11;  %.5;   % mm    %WORKS BRILLIANTLY

lens = lensC('apertureSample', [nSamples nSamples], ...
    'fileName', lensFileName, ...
    'apertureMiddleD', apertureMiddleD, ...
    'diffractionEnabled', true);

% Create a film (sensor) 
% position - relative to center of final lens surface
% size - 'mm'
% wavelength samples
lens.set('wave', [400:10:700 ]);

wave = lens.get('wave');

%put it 16 mm away
film = filmC('position', [0 0 50], ...
    'resolution', [300 300 1], ...
    'size', [2/sqrt(2) 2/sqrt(2)], ...
    'wave', wave);

% Create a camera out of lens, film ,and point source
camera = psfCameraC('lens',lens,'film',film,'point source',point{1});

% Sequence of events for estimating the PSF, 
nLines = 100;
jitter = false;
%camera.estimatePSF(nLines,jitter);

%limits the entrance aperture so this can run faster
subsection = [];

%method = 'HURB';
method = 'huygens';
rtType = 'ideal';
camera.estimatePSF(nLines,jitter,subsection, method, rtType);

oiHuygens = camera.oiCreate(); vcAddObject(oiHuygens); oiWindow;
    

%% Produce Theoretical results
%scene = sceneCreate('point array',256,128);
%scene = sceneSet(scene,'fov',8);
%vcAddAndSelectObject(scene); sceneWindow;

%load scene file
d = displayCreate('equal energy');
scene = sceneFromFile('pointTest.png', 'rgb', [], d);
wave = sceneGet(scene, 'wave');
onesPhotons = ones(size(wave)) * 10^15;
equalPhotonsEnergy = Quanta2Energy(wave, onesPhotons);
scene = sceneAdjustIlluminant(scene, equalPhotonsEnergy);

%assign parameters from above
sensorWidth = film.size(1);
filmDistance = film.position(3);
focalLength = 50;
apertureDiameter = apertureMiddleD;

horFieldofView = 2 * atan(sensorWidth/(2 * filmDistance)) * 180/pi * .8;
scene = sceneSet(scene,'fov',horFieldofView);
scene = sceneSet(scene, 'distance', 100000001);   %scene = sceneSet(scene, 'distance', 2001);
vcAddAndSelectObject(scene); sceneWindow;

%create optical image
oiT = oiCreate;
optics = oiGet(oiT,'optics'); 
fNumber = focalLength/apertureDiameter;
optics = opticsSet(optics,'fnumber',fNumber);

% In this example we set the properties of the optics to include cos4th
% falloff for the off axis vignetting of the imaging lens
optics = opticsSet(optics,'offaxis','cos4th');
optics = opticsSet(optics,'focallength',focalLength * 10^-3);    
oiT = oiSet(oiT,'optics',optics);
oiT = oiCompute(scene,oiT);
vcAddAndSelectObject(oiT); oiWindow;

%% Plot mesh plots for all 3 techniques

%Theoretical
oiPhotons = oiGet(oiT, 'photons');
PSFLineSpectral = sum(oiPhotons, 1);
PSFLineSpectral = reshape(PSFLineSpectral, [size(oiPhotons,1) size(oiPhotons, 3)]);
plotBound = sensorWidth/2 * 10^3;

[X, Y] = meshgrid(400:10:700, linspace(plotBound, -plotBound, size(PSFLineSpectral, 1)));
figure; mesh(X, Y, PSFLineSpectral./max(PSFLineSpectral(:)));
xlabel('Wavelength (nm)')
ylabel('Position (um)')
zlabel('Intensity (rel.)');
title('Theoretical Linespread');


%HURB
oiPhotons = oiGet(oiHURB, 'photons');
PSFLineSpectral = sum(oiPhotons, 1);
PSFLineSpectral = reshape(PSFLineSpectral, [size(oiPhotons,1) size(oiPhotons, 3)]);
plotBound = sensorWidth/2 * 10^3;

[X, Y] = meshgrid(400:10:700, linspace(plotBound, -plotBound, size(PSFLineSpectral, 1)));
figure; mesh(X, Y, PSFLineSpectral./max(PSFLineSpectral(:)));
xlabel('Wavelength (nm)')
ylabel('Position (um)')
zlabel('Intensity (rel.)');
title('HURB Linespread');

%Huygens
oiPhotons = oiGet(oiHuygens, 'photons');
PSFLineSpectral = sum(oiPhotons, 1);
PSFLineSpectral = reshape(PSFLineSpectral, [size(oiPhotons,1) size(oiPhotons, 3)]);
plotBound = sensorWidth/2 * 10^3;

[X, Y] = meshgrid(400:10:700, linspace(plotBound, -plotBound, size(PSFLineSpectral, 1)));
figure; mesh(X, Y, PSFLineSpectral./max(PSFLineSpectral(:)));
xlabel('Wavelength (nm)')
ylabel('Position (um)')
zlabel('Intensity (rel.)');  
title('Huygens-Fresnel Linespread');

%% plot line vs. wavelength plots (linespread) plot 3 PSFs on 1 figure
oiPhotonsTemp = oiGet(oiT, 'photons');
PSFLineT = sum( oiPhotonsTemp(:,:,16), 1);
PSFLineTS = PSFLineT /max(PSFLineT(:));

oiPhotonsTemp = oiGet(oiHuygens, 'photons');
PSFLineHuygens = sum(oiPhotonsTemp(:,:, 16) , 1);
PSFLineHuygens = PSFLineHuygens / max(PSFLineHuygens(:));

% oiPhotonsTemp = oiGet(oiHURBTuned, 'photons');
% PSFLineHURBTuned = sum(oiPhotonsTemp(:,:,16), 1);
% PSFLineHURBTuned = PSFLineHURBTuned / max(PSFLineHURBTuned);

oiPhotonsTemp = oiGet(oiHURB, 'photons');
PSFLineHURB = sum(oiPhotonsTemp(:,:,16), 1);
PSFLineHURB = PSFLineHURB / max(PSFLineHURB);

positionT = linspace(-sensorWidth/2 *1000, sensorWidth/2 *1000, length(PSFLineT));
position = linspace(-sensorWidth/2 * 1000, sensorWidth/2 * 1000, length(PSFLineHURB));
figure;
plot( positionT, PSFLineTS, position, PSFLineHURB, position, PSFLineHuygens);


title(['Linespread Comparison at 550nm;' num2str(focalLength) 'mm;f/' ...
    num2str(focalLength/apertureDiameter, 2)  ]);
xlabel('um')
%axis([-40 40 0 1]);  %don't show the bad part of the theoretical plot
ylabel('Relative radiance');
legend('Theoretical', 'HURB', 'Huygens-Fresnel');
% 
% 
% %save figure as a tiff file
% fileName = ['PSFC_' num2str(sampleArray{index}.focalLength) 'mm_f' ...
%     num2str(focalLength/(apertureDiameter))];
% hgexport(gcf, [fileName '.tif'], hgexport('factorystyle'), 'Format', 'tiff');

%% plot line vs. wavelength plots (PSF slice) plot 3 PSFs on 1 figure
oiPhotonsTemp = oiGet(oiT, 'illuminance');
PSFLineT = oiPhotonsTemp((size(oiPhotonsTemp, 1))/2,:);
PSFLineTS = PSFLineT /max(PSFLineT(:));

oiPhotonsTemp = oiGet(oiHuygens, 'illuminance');
PSFLineHuygens = oiPhotonsTemp((size(oiPhotonsTemp, 1))/2,:);
PSFLineHuygens = PSFLineHuygens / max(PSFLineHuygens(:));

oiPhotonsTemp = oiGet(oiHURB, 'illuminance');
PSFLineHURB = oiPhotonsTemp((size(oiPhotonsTemp, 1))/2,:);
PSFLineHURB = PSFLineHURB / max(PSFLineHURB);

positionT = linspace(-sensorWidth/2 *1000, sensorWidth/2 *1000, length(PSFLineT));
position = linspace(-sensorWidth/2 * 1000, sensorWidth/2 * 1000, length(PSFLineHURB));
figure;
plot( positionT, PSFLineTS, position, PSFLineHURB, position, PSFLineHuygens);


title(['PSF Slice Comparison at 550nm;' num2str(focalLength) 'mm;f/' ...
    num2str(focalLength/apertureDiameter, 2)  ]);
xlabel('um')
%axis([-40 40 0 1]);  %don't show the bad part of the theoretical plot
ylabel('Relative radiance');
legend('Theoretical', 'HURB', 'Huygens-Fresnel');