%% s_s3dRenderDesk.m
%
% This should render the desk scene for flash/noflash cases.
%
% We should use the basic format in bench and metronome files to do the
% rendering, but for the desk/default.pbrt.

%%
ieInit

%% PBRT will run the PBRT script
pbrtDir = fullfile(s3dRootPath, 'data', 'pbrtScenes','desk');
if ~exist(pbrtDir,'dir'), error('Desk pbrt directory is missing'); end

% Governs the number of samples per pixel. A minimum is 32.
% For the scarlet one we ran it with 4096.s
nSamples = 1024; %512;

%%

%initialization
clear curPbrt;
curPbrt = pbrtObject();

% Specify properties files (materials and geometry)
matFile = fullfile(pbrtDir, 'default-mat.pbrt');
geoFile = fullfile(pbrtDir, 'default-geom.pbrt');


%% Standard camera properties

% This is from the default.pbrt file for desk.
LookAt  = ...
[ -127.5649  -89.6985   80.9301;
   -81.3417  -63.8818   72.3964;
     0 0 1];
curPbrt.camera.setPosition(LookAt);

% light properties
spectrum = pbrtSpectrumObject('rgb I', [1000 1000 1000]);
% lightFrom = [  -56.914787 -105.385544 35.0148];  % Position of source
% lightTo =   [-56.487434 -104.481461 34.8  ];       % Direction of principal ray
coneAngle      = 180;    % Angle of rays from light source
coneDeltaAngle = 180;    % Drop off of illumination???

lightSource = pbrtLightSpotObject('light', spectrum, coneAngle, coneDeltaAngle, LookAt(1,:), LookAt(2,:));  


%% specify the wide field camera 
    
% From the default.pbrt file
curPbrt.camera.setLens(pbrtLensRealisticObject());   
curPbrt.camera.lens.filmDistance = 52; % 70;  % 133.33;
curPbrt.camera.lens.filmDiag = 43.3;
curPbrt.camera.lens.specFile = 'dgauss.50mm.dat';
curPbrt.camera.lens.apertureDiameter = 3; % in mm
curPbrt.camera.lens.curveRadius = 0;       % Experimental
curPbrt.camera.setResolution(450, 300);
% curPbrt.camera.setResolution(64, 32);

%% add old parts, put in new ones
% Candidate for a function
curPbrt.removeMaterial();
curPbrt.addMaterial(matFile);
curPbrt.removeGeometry();
curPbrt.addGeometry(geoFile);
curPbrt.removeLight();
curPbrt.addLightSource(lightSource);

% set sampler
curPbrt.sampler.removeProperty();
curPbrt.sampler.addProperty(pbrtPropertyObject('integer pixelsamples', nSamples));

%% Render the oi
oiName = 'desk';
dockerFlag = true;
oi = s3dRenderOIAndDepthMap(curPbrt, oiName, dockerFlag);

vcAddAndSelectObject(oi);
oiWindow;

%% Deprecated.
% We should use the modern methods to make this scene.
% 
% 
% mkdir('tempOutput');
% chdir('tempOutput');
% unix('rm *');
% 
% %% scene rendering
% % unix([fullfile(pbrtHome, '/src/bin/pbrt') fname '--outfile output.dat']);
% outfile = 'desk_out.dat';
% unix([fullfile(pbrtHome, '/src/bin/pbrt ') fname ' --outfile ' outfile]);
% 
% % ISET will read the PBRT output
% % scene = sceneSet(scene,'fov', 8);
% oi = pbrt2oi(outfile);
% % oi = oiSet (oi, 'horizontalfieldofview', 8 * 200/150 );
% vcAddAndSelectObject(oi);
% oiWindow;
% m = oiGet(oi, 'mean illuminance')
% 
% %% camera processing - no flash image
% load('deskNoFlashOiRight.mat');
% oi = opticalimage;
% oi = oiSet(oi, 'photons', oiGet(oi,'photons') * .5 * 10^14);  %some normalization issues
% myOptics = oiGet(oi, 'optics');  %to create proper crop at sensor
% myOptics = opticsSet(myOptics,'focalLength', .050);  %distance from lens to sensor - sensor will inherit this value!
% oi = oiSet(oi, 'optics', myOptics);
% oi = oiSet(oi,'fov', 39.6); %fov of the scene
% % we want this to MATCH the fov of the sensor, so that no cropping occurs
% vcAddAndSelectObject(oi); oiWindow;
% 
% %sensor processing
% %sensor = s3dProcessSensor(oi, .0096, [600 400]);  
% % sensor = s3dProcessSensor(oi, .01, [600 400], .04);    %why does exposure need to be set to this level?  what does this mean?  is pink overexposure realistic?
% sensor = s3dProcessSensor(oi, [], [600 400], .04);    %why does exposure need to be set to this level?  what does this mean?  is pink overexposure realistic?
% vcAddAndSelectObject('sensor',sensor); sensorImageWindow;
% 
% %image processing
% image = s3dProcessImage(sensor);
% vcAddAndSelectObject(image); vcimageWindow;
% 
% 
% %% camera processing -  flash image
% load('deskFlashOi.mat');
% oi = opticalimage;
% oi = oiSet(oi, 'photons', oiGet(oi,'photons') * 10^9);  %some normalization issues
% myOptics = oiGet(oi, 'optics');  %to create proper crop at sensor
% myOptics = opticsSet(myOptics,'focalLength', .050);  %distance from lens to sensor - sensor will inherit this value!
% oi = oiSet(oi, 'optics', myOptics);
% oi = oiSet(oi,'fov', 39.6); %fov of the scene   %26.99
% % we want this to MATCH the fov of the sensor, so that no cropping occurs
% vcAddAndSelectObject(oi); oiWindow;
% 
% %sensor processing
% sensor = s3dProcessSensor(oi, [], [600 400], .010);  
% vcAddAndSelectObject('sensor',sensor); sensorImageWindow;
% 
% %image processing
% image = s3dProcessImage(sensor);
% vcAddAndSelectObject(image); vcimageWindow;
% 





%% - old stuff

% 
% 
% 
% %initialize
% outfile = 'desk_out.dat';   %this name may be changed for clarity purposes
% chdir(fullfile(s3dRootPath, 'scripts', 'pbrtFiles'));
% 
% %call pbrt unix command
% unix([fullfile(pbrtHome, '/src/bin/pbrt') ' desk/default.pbrt --outfile ' outfile]);
% 
% % ISET will read the PBRT output
% % scene = sceneSet(scene,'fov', 8);
% oi = pbrt2oi(outfile);
% oi = oiSet(oi, 'photons', oiGet(oi,'photons') * 10^14);  %some normalization issues
% myOptics = oiGet(oi, 'optics');  %to create proper crop at sensor
% myOptics = opticsSet(myOptics,'focalLength', .050);  %distance from lens to sensor - sensor will inherit this value!
% oi = oiSet(oi, 'optics', myOptics);
% oi = oiSet(oi,'fov', 38.07); %fov of the scene
% % we want this to MATCH the fov of the sensor, so that no cropping occurs
% vcAddAndSelectObject(oi); oiWindow;


