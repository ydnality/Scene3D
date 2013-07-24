% Runs PBRT and imports it in ISET for the bench scene. 

%% PBRT will run the PBRT script - initializing
chdir(fullfile(s3dRootPath, 'scripts', 'pbrtFiles'));
% fname = '../indestructibleObject/default.pbrt'
fname = '../indestructibleObject/lambertian-down.pbrt'
fname = '../indestructibleObject/textured-down.pbrt'
fname = '../indestructibleObject/graycard-down.pbrt'
fname = '../indestructibleObject/simpleReflectance-down.pbrt'

mkdir('tempOutput');
chdir('tempOutput');
unix('rm *');

% scene rendering
% unix([fullfile(pbrtHome, '/src/bin/pbrt') fname '--outfile output.dat']);
outfile = 'indestructibleObject_out.dat';
dMapFile = 'indestructibleObject_out_DM.dat'; 
unix([fullfile(pbrtHome, '/src/bin/pbrt ') fname ' --outfile ' outfile]);

% ISET will read the PBRT output
% scene = sceneSet(scene,'fov', 8);
oi = pbrt2oi(outfile);
% oi = oiSet (oi, 'horizontalfieldofview', 8 * 200/150 );
vcAddAndSelectObject(oi);
oiWindow;
m = oiGet(oi, 'mean illuminance')


% Runs PBRT and imports it in ISET for the bench scene. 


%%  read and output depth map
% ** make sure the rendering file has a small initial aperture, and only 1
% sample per pixel!!!

outfile = 'indObj_out.dat';
dMapFile = 'indObj_out_DM.dat';
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

% this code is only needed when doing multiple samples... we don't need
% that for a realistic camera
% numSamples = 1;
% ratio = 450/300;
% for j = 1:size(depthMap,1)
%     for i = round(j*ratio):size(depthMap,2)
%         depthMap(j,i) = depthMap(j,i)/numSamples;
%     end
% end


% figure; imagesc(depthMap);



%% TEST CASE - Walls straight on ... PBRT will run the PBRT script
chdir(fullfile(s3dRootPath, 'scripts', 'pbrtFiles'));
% fname = '../indestructibleObject/default.pbrt'
fname = '../floorWallBottomBack/lambertian.pbrt'

mkdir('tempOutput');
chdir('tempOutput');
unix('rm *');

%scene rendering
% unix([fullfile(pbrtHome, '/src/bin/pbrt') fname '--outfile output.dat']);
outfile = 'floorWallBottomBack_out.dat';
% dMapFile = 'indestructibleObject_out_DM.dat'; 
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

outfile = 'floorWallBottomBack_out.dat';
dMapFile = 'floorWallBottomBack_out_DM.dat';

outfile = 'indObj_out.dat';
dMapFile = 'indObj_out_DM.dat';
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

% this code is only needed when doing multiple samples... we don't need
% that for a realistic camera
% numSamples = 1;
% ratio = 450/300;
% for j = 1:size(depthMap,1)
%     for i = round(j*ratio):size(depthMap,2)
%         depthMap(j,i) = depthMap(j,i)/numSamples;
%     end
% end


% figure; imagesc(depthMap);



%% camera processing - no flash image
load('indObjNoFlashOi.mat'); %load('indObjNoFlashOi.mat');

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
% load('indObjFlashOiLambertianNoAmbient.mat');
% load('indObject2FlashDepth/frontFlashOiLambertian.mat');
% load('indObject2FlashDepth/grayBackOi.mat');  % currently testing
load('indObject2FlashDepth/downBackFlashOi.mat');  % currently testing

% load('floorWallBottomBack/backFlashDown100Oi.mat');
% load('floorWallBottomBack/sideFlashDownOi.mat');
% load('floorWallBottomBack/sideFlashDown25Oi.mat');

% load('indObject2FlashDepth/backFlashOiLambertianCloser.mat');
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
% sensor = s3dProcessSensor(oi, 0, [], .03, []);      %front auto-exposure time:
%autoExpTime = .4768
%sensor = s3dProcessSensor(oi, 0, [], [], .4768); 
% sensor = s3dProcessSensor(oi, 0, [], [], .05);   %ind object
% sensor = s3dProcessSensor(oi, 0, [], [], .04);    %test scene - front
% flash
% sensor = s3dProcessSensor(oi, 0, [], [], .08);    %test scene - back flash - multiplication factor - 2


% sensor = s3dProcessSensor(oi, 0, [], [], .64);   %test scene - back 100m
% flash - multiplication factor - 16 
sensor = s3dProcessSensor(oi, 0, [],  .32);   %test scene - back 100m flash - multiplication factor - 8
% sensor = s3dProcessSensor(oi, 0, [], [], .02);    %multiplication factor - .5

% sensor = s3dProcessSensor(oi, 0, [], [], .05); 
% sensor = s3dProcessSensor(oi, 0, [], [], .002); 
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

%image processing
% [image, transformMatrices] = s3dProcessImage(sensor, []);  %first time running


% wantedTransformMatrix = ;  % need to find this first... 
[image, transformMatrices] = s3dProcessImage(sensor, []);

vcAddAndSelectObject(image); vcimageWindow;




%% new depth-map method for perspective camera - consider putting this into a function

% [path,name,ext] = fileparts(outfile);
% dMapFile = [path '/' name '_DM.dat'];
% badDepthMap = s3dReadDepthMapFile(dMapFile);
% newDepthMap = badDepthMap;
% numSamples = 64;
% 
% for j = 1:size(badDepthMap,1);
%     for i = j+1:size(badDepthMap,2);
%         newDepthMap(j, i) =  badDepthMap(j, i)/numSamples;
%     end
% end
% 
% figure; imagesc(badDepthMap);
% figure; imagesc(newDepthMap);

%% depth map rendering experiment - old version - do not run

%basically, this renders a set number of depth maps, then combines them
%together.  They are combined using something that is NOT an averaging
%technique, such as using the median operation, or minimum operation.
%Initial evaluation shows that the median method works best.

%We plan to use rendertoolbox to make this rendering technique more
%streamlined and elegant in the future.

% % % numDM =  101;  %number of depth map files to render
% % % depthMap = zeros(100, 100, numDM);
% % % for i = 1:numDM    
% % %     % depth map rendering
% % %     [path,name,ext] = fileparts(fname); 
% % %    
% % %     unix([fullfile(pbrtHome, '/src/bin/pbrt ') path '/' name '_DM.pbrt --outfile depth.dat' ]);
% % %     
% % %     dMapFile = 'depth_DM.dat'; 
% % %     depthMap(:,:,i) = s3dReadDepthMapFile(dMapFile);
% % %     unix('rm *');
% % % end
% % % 
% % % depthMapProcessedMin = min(depthMap, [], 3);
% % % depthMapProcessedMedian = median(depthMap, 3);
% % % figure; imagesc(depthMap(:,:,1));
% % % figure; imagesc(depthMapProcessedMin);
% % % figure; imagesc(depthMapProcessedMedian);
% % % 
% % % 
% % % %display the depth map
% % % %figure; imagesc(oi.depthMap);
% % % unix('cd ..');


% 
% %% ISET will read the PBRT output
% % scene = sceneSet(scene,'fov', 8);
% %oi = pbrt2oi('benchScene_dca_zoom_hq.dat');
% % oi = pbrt2oi('benchScene_d_zoom_hq.dat');
% oi = pbrt2oi('benchScene_nod_zoom_hq.dat');
% 
% 
% % oi = oiSet (oi, 'horizontalfieldofview', 8 * 200/150 );
% vcAddAndSelectObject(oi);
% oiWindow;
% 
% m = oiGet(oi, 'mean illuminance')
% unix('cd ..');



