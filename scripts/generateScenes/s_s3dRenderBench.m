

% Runs PBRT and imports it in ISET for the bench scene. 

%% PBRT will run the PBRT script
chdir(fullfile(s3dRootPath, 'scripts', 'pbrtFiles'));

% rendering for diffraction/chromatic aberrations
unix([fullfile(pbrtHome, '/src/bin/pbrt') ' benchScene/defaultBiggerZoom.pbrt']);

%% ISET will read the PBRT output
% scene = sceneSet(scene,'fov', 8);
oi = pbrt2oi('benchScene_d.dat');
% oi = oiSet (oi, 'horizontalfieldofview', 8 * 200/150 );
vcAddAndSelectObject(oi);
oiWindow;

m = oiGet(oi, 'mean illuminance')
unix('cd ..');




%% no-flash rendering

unix([fullfile(pbrtHome, '/src/bin/pbrt') ' benchScene/defaultBiggerZoom_NoFlash.pbrt']);
% ISET will read the PBRT output
% scene = sceneSet(scene,'fov', 8);
oi = pbrt2oi('benchScene_d.dat');
% oi = oiSet (oi, 'horizontalfieldofview', 8 * 200/150 );
vcAddAndSelectObject(oi);
oiWindow;
m = oiGet(oi, 'mean illuminance')


%% flash rendering
unix([fullfile(pbrtHome, '/src/bin/pbrt') ' benchScene/defaultBiggerZoom_Flash.pbrt']);
% ISET will read the PBRT output
% scene = sceneSet(scene,'fov', 8);
oi = pbrt2oi('benchScene_d.dat');
% oi = oiSet (oi, 'horizontalfieldofview', 8 * 200/150 );
vcAddAndSelectObject(oi);
oiWindow;
m = oiGet(oi, 'mean illuminance');
unix([fullfile(pbrtHome, '/src/bin/pbrt') ' benchScene/defaultBiggerZoom_NoFlash.pbrt --outfile benchScene_d_noflash.dat']);

% ISET will read the PBRT output
% scene = sceneSet(scene,'fov', 8);
oi = pbrt2oi('benchScene_d_noflash.dat');
oi = oiSet(oi, 'photons', oiGet(oi,'photons') * 10^14);  %some normalization issues
myOptics = oiGet(oi, 'optics');  %to create proper crop at sensor
myOptics = opticsSet(myOptics,'focalLength', .050);  %distance from lens to sensor - sensor will inherit this value!
oi = oiSet(oi, 'optics', myOptics);
oi = oiSet(oi,'fov', 38.07); %fov of the scene
% we want this to MATCH the fov of the sensor, so that no cropping occurs
vcAddAndSelectObject(oi); oiWindow;

%sensor processing
sensor = s3dProcessSensor(oi, .0096);  
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

%image processing
image = s3dProcessImage(sensor);
vcAddAndSelectObject(image); vcimageWindow;














%% flash rendering
unix([fullfile(pbrtHome, '/src/bin/pbrt') ' benchScene/defaultBiggerZoom_Flash.pbrt --outfile benchScene_d_flash.dat']);
oi = pbrt2oi('benchScene_d_flash.dat');
oi = oiSet(oi, 'photons', oiGet(oi,'photons') * 10^14);  %some normalization issues
myOptics = oiGet(oi, 'optics');  %to create proper crop at sensor
myOptics = opticsSet(myOptics,'focalLength', .050);  %distance from lens to sensor - sensor will inherit this value!
oi = oiSet(oi, 'optics', myOptics);
oi = oiSet(oi,'fov', 39.6); %fov of the scene
% we want this to MATCH the fov of the sensor, so that no cropping occurs
vcAddAndSelectObject(oi); oiWindow;, 

%sensor processing
sensor = s3dProcessSensor(oi, .000096, [400 600], .025);  
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

%image processing
image = s3dProcessImage(sensor);
vcAddAndSelectObject(image); vcimageWindow;


