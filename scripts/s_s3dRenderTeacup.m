% renders and displays the desk scene

%% PBRT will run the PBRT script
chdir(fullfile(s3dRootPath, 'scripts', 'pbrtFiles'));
fname = '../teacup/teacup.pbrt'

mkdir('tempOutput');
chdir('tempOutput');
unix('rm *');

%% scene rendering
% unix([fullfile(pbrtHome, '/src/bin/pbrt') fname '--outfile output.dat']);
outfile = 'teacup_out.dat';
unix([fullfile(pbrtHome, '/src/bin/pbrt ') fname ' --outfile ' outfile]);

% ISET will read the PBRT output
% scene = sceneSet(scene,'fov', 8);
oi = pbrt2oi(outfile);
% oi = oiSet (oi, 'horizontalfieldofview', 8 * 200/150 );
vcAddAndSelectObject(oi);
oiWindow;
m = oiGet(oi, 'mean illuminance')

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
% 