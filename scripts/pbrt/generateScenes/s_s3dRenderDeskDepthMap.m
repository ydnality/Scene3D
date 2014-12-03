% Runs PBRT and imports it in ISET for the bench scene. 

%% PBRT will run the PBRT script
chdir(fullfile(s3dRootPath, 'scripts', 'pbrtFiles'));
fname = '../desk/default.pbrt'

mkdir('tempOutput');
chdir('tempOutput');
unix('rm *');

%% scene rendering
% unix([fullfile(pbrtHome, '/src/bin/pbrt') fname '--outfile output.dat']);
outfile = 'desk_out.dat';
unix([fullfile(pbrtHome, '/src/bin/pbrt ') fname ' --outfile ' outfile]);

% ISET will read the PBRT output
% scene = sceneSet(scene,'fov', 8);
oi = pbrt2oi(outfile);
% oi = oiSet (oi, 'horizontalfieldofview', 8 * 200/150 );
vcAddAndSelectObject(oi);
oiWindow;
m = oiGet(oi, 'mean illuminance')

%% new depth-map method for perspective camera - consider putting this into a function

[path,name,ext] = fileparts(outfile);
dMapFile = [path '/' name '_DM.dat'];
badDepthMap = s3dReadDepthMapFile(dMapFile);
newDepthMap = badDepthMap;
numSamples = 128;

for j = 1:size(badDepthMap,1);
    for i = j+1:size(badDepthMap,2);
        newDepthMap(j, i) =  badDepthMap(j, i)/numSamples;
    end
end

figure; imagesc(badDepthMap);
figure; imagesc(newDepthMap);

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



