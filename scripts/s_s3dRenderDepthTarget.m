%% PBRT will run the PBRT script - initializing
chdir(fullfile(s3dRootPath, 'scripts', 'pbrtFiles'));
fname = '../depthTarget/depthTarget.pbrt'
fname = '../depthTarget/depthTargetHorizontal.pbrt'
fname = '../depthTargetGradient/depthTargetGradient.pbrt'
fname = '../depthTargetDepths/depthTargetDepths.pbrt'

mkdir('tempOutput');
chdir('tempOutput');
unix('rm *');

% scene rendering
outfile = 'depthTarget_out.dat';
dMapFile = 'depthTarget_out_DM.dat'; 
unix([fullfile(pbrtHome, '/src/bin/pbrt ') fname ' --outfile ' outfile]);

% ISET will read the PBRT output
% scene = sceneSet(scene,'fov', 8);
oi = pbrt2oi(outfile);
% oi = oiSet (oi, 'horizontalfieldofview', 8 * 200/150 );
vcAddAndSelectObject(oi);
oiWindow;
m = oiGet(oi, 'mean illuminance')


%% render depth map
% read and output depth map
% ** make sure the rendering file has a small initial aperture, and only 1
% sample per pixel!!!

outfile = 'depthTarget_out.dat';
dMapFile = 'depthTarget_out_DM.dat'; 

imageWidth = 450;
imageHeight = 300;
numRenders = 31;
depthMap = zeros(imageHeight, imageWidth, numRenders);

for i = 1:numRenders
    
    unix([fullfile(pbrtHome, '/src/bin/pbrt ') fname ' --outfile ' outfile]);
    depthMap(:,:, i) = s3dReadDepthMapFile(dMapFile, [300 450]);
    unix('rm *');
end

depthMapProcessedMedian = median(depthMap, 3);
figure; imagesc(depthMapProcessedMedian);


%% camera processing - no flash image
load('depthTarget/noFlashOi.mat'); 

oi = opticalimage;
oi = oiSet(oi, 'photons', oiGet(oi,'photons') * 10^13);  %some normalization issues
myOptics = oiGet(oi, 'optics');  %to create proper crop at sensor
myOptics = opticsSet(myOptics,'focalLength', .050);  %distance from lens to sensor - sensor will inherit this value!
oi = oiSet(oi, 'optics', myOptics);
oi = oiSet(oi,'fov', 39.60); %fov of the scene
% we want this to MATCH the fov of the sensor, so that no cropping occurs
vcAddAndSelectObject(oi); oiWindow;

%sensor processing
% sensor = s3dProcessSensor(oi, .0096, [], .03);     %high noise
sensor = s3dProcessSensor(oi, 0, [], 0);    %low noise
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

%image processing
image = s3dProcessImage(sensor);
vcAddAndSelectObject(image); vcimageWindow;


%% camera processing -  flash image
load('depthTarget/flashOi.mat')  

oi = opticalimage;
% oi = oiSet(oi, 'photons', oiGet(oi,'photons') * 10^14);  %some normalization issues
oi = oiSet(oi, 'photons', oiGet(oi,'photons') * 10^13);  %some normalization issues
myOptics = oiGet(oi, 'optics');  %to create proper crop at sensor
myOptics = opticsSet(myOptics,'focalLength', .050);  %distance from lens to sensor - sensor will inherit this value!
oi = oiSet(oi, 'optics', myOptics);
oi = oiSet(oi,'fov', 39.60); %fov of the scene
% we want this to MATCH the fov of the sensor, so that no cropping occurs
vcAddAndSelectObject(oi); oiWindow;

%sensor processing
sensor = s3dProcessSensor(oi, 0, [], [], .32);   %test scene - back 100m flash - multiplication factor - 8
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

%image processing
[image, transformMatrices] = s3dProcessImage(sensor, []);
vcAddAndSelectObject(image); vcimageWindow;

