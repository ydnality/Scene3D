% Use the point test file, blur it in 2 different ways.  The first way is
% to use the new ray tracing method which uses Heisenburg Uncertainty Ray
% Bending (HURB).  The second way is the classical way, using theoretical PSF's.  

%% HURB ray tracing results

% s_s3dRenderDiffractionPoint;  %use this if no prerendered images to load

%for this example only
% load 'pointTestOI.mat';
% oi = opticalimage;
%load 'generatedDiffractioNoSpectral.mat'
load 'generatedDiffraction101612.mat'
oi = opticalimage;

oiIlluminance = oiGet(oi, 'illuminance');
PSFLine = oiIlluminance(size(oiIlluminance,1)/2, :);
oi = oiSet (oi, 'horizontalfieldofview', 8 * 200/150 );  %need to correct for "horizontalness"
vcAddAndSelectObject(oi);

position = linspace(-359.5, 359.5, length(PSFLine));
figure;
plot(position, PSFLine);
title('Raytracing PSF');
xlabel('um');
ylabel('Illuminance');


%% Theoretical results

%scene = sceneCreate('point array',256,128);
%scene = sceneSet(scene,'fov',8);

%vcAddAndSelectObject(scene); sceneWindow;

%load scene file
 scene = sceneFromFile('pointTest.png', 'rgb');
%scene = sceneFromFile('usairforce.png', 'rgb');

scene = sceneSet(scene,'fov',8);
scene = sceneSet(scene, 'distance', 2);
vcAddAndSelectObject(scene); sceneWindow;

%create optical image
oiT = oiCreate;
optics = oiGet(oiT,'optics');           
optics = opticsSet(optics,'fnumber',22);
% In this example we set the properties of the optics to include cos4th
% falloff for the off axis vignetting of the imaging lens
optics = opticsSet(optics,'offaxis','cos4th');
optics = opticsSet(optics,'focallength',2.2e-3);    
oiT = oiSet(oiT,'optics',optics);
oiT = oiCompute(scene,oiT);
vcAddAndSelectObject(oiT); oiWindow;

%plot PSF
oiIlluminanceT = oiGet(oiT, 'illuminance');
PSFLineT = oiIlluminanceT(size(oiIlluminanceT,1)/2, :);
PSFLineTS = PSFLineT * max(PSFLine(:))/max(PSFLineT(:));
positionT = linspace(-192.2, 192.2, length(PSFLineT));
figure;
plot(position, PSFLine, positionT, PSFLineTS);
title('Theoretical PSF');
xlabel('um');
ylabel('Illuminance');