% Runs PBRT and imports it in ISET for the bench scene. 

%% PBRT will run the PBRT script
chdir(fullfile(s3dRootPath, 'scripts', 'pbrtFiles'));
% fname = '../indestructibleObject/default.pbrt'
fname = '../indestructibleObject/lambertian.pbrt'

mkdir('tempOutput');
chdir('tempOutput');
unix('rm *');

%% scene rendering
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

%% render depth map
%read and output depth map
outfile = 'indestructibleObject_out.dat';
dMapFile = 'indestructibleObject_out_DM.dat'; 

unix([fullfile(pbrtHome, '/src/bin/pbrt ') fname ' --outfile ' outfile]);
depthMap = s3dReadDepthMapFile(dMapFile, [300 450]);
numSamples = 4096;

ratio = 450/300;
for j = 1:size(depthMap,1)
    for i = round(j*ratio):size(depthMap,2)
        depthMap(j,i) = depthMap(j,i)/numSamples;
    end
end
figure; imagesc(depthMap);



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
load('indObject2FlashDepth/grayBackOi.mat');
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
sensor = s3dProcessSensor(oi, 0, [], [], .05); 
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



