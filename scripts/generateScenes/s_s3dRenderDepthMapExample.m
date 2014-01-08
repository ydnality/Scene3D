% Runs PBRT and imports it in ISET for the bench scene. 

%% PBRT will run the PBRT script
chdir(fullfile(s3dRootPath, 'scripts', 'pbrtFiles'));
fname = ' ../cones/defaultBiggerZoom.pbrt'

mkdir('tempOutput');
chdir('tempOutput');
unix('rm *');

%% cones scene rendering
unix([fullfile(pbrtHome, '/src/bin/pbrt') fname]);

% ISET will read the PBRT output
% scene = sceneSet(scene,'fov', 8);
oi = pbrt2oi('output_d.dat');
% oi = oiSet (oi, 'horizontalfieldofview', 8 * 200/150 );
vcAddAndSelectObject(oi);
oiWindow;
m = oiGet(oi, 'mean illuminance')

%% depth map rendering experiment

%basically, this renders a set number of depth maps, then combines them
%together.  They are combined using something that is NOT an averaging
%technique, such as using the median operation, or minimum operation.
%Initial evaluation shows that the median method works best.

%We plan to use rendertoolbox to make this rendering technique more
%streamlined and elegant in the future.

numDM =  31;  %number of depth map files to render
depthMap = zeros(300, 300, numDM);
for i = 1:numDM    
    % depth map rendering
    [path,name,ext] = fileparts(fname); 
   
    unix([fullfile(pbrtHome, '/src/bin/pbrt') path '/' name '_DM.pbrt' ]);
    
    dMapFile = 'output_d_dm_DM.dat'; 
    depthMap(:,:,i) = s3dReadDepthMapFile(dMapFile);
    unix('rm *');
end

depthMapProcessedMin = min(depthMap, [], 3);
depthMapProcessedMedian = median(depthMap, 3);
figure; imagesc(depthMap(:,:,1));
figure; imagesc(depthMapProcessedMin);
figure; imagesc(depthMapProcessedMedian);


%display the depth map
%figure; imagesc(oi.depthMap);
unix('cd ..');


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
