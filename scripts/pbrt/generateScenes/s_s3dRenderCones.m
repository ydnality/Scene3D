% Runs PBRT and imports it in ISET for the cones scene. 

%% This section runs the PBRT script directly (demo) 
dockerFlag = false;
fname = fullfile(dataPath, 'pbrtScenes', 'cones', 'defaultBiggerZoom.pbrt');  %defaultBiggerZoom_depth.pbrt is for depth map
oi = s3dRenderOIAndDepthMap(fname, [], 'cones', dockerFlag);
vcAddObject(oi); oiWindow;

%% This section uses a PBRT wrapper for more flexibility for rendering

%initialization
clear curPbrt;
curPbrt = pbrtObject();

includeFile = fullfile(dataPath, 'pbrtScenes', 'cones', 'cones.pbrt');
curPbrt.addInclude(includeFile);
curPbrt.camera.addTransform(pbrtTransformObject('Scale', [1000 1000 1000]));
curPbrt.camera.addTransform(pbrtTransformObject('Rotate',  [-3 1 0 0]));
curPbrt.camera.addTransform(pbrtTransformObject('Rotate', [52 0 1 0]));
curPbrt.camera.addTransform(pbrtTransformObject('Translate', [-2.3 -.05 .5]));

%make a lens and add to pbrt object
filmDist = 36.77;
filmDiag = 10;
specFile = 'dgauss.50mm.dat';
apertureDiameter = 3;
diffraction = false;
chromaticAberration =false;
lens = pbrtLensRealisticObject(filmDist, filmDiag, specFile, apertureDiameter, diffraction, chromaticAberration, []);
curPbrt.camera.setLens(lens);

%make an area light and add to pbrt object
light = pbrtAreaLightObject('area', pbrtSpectrumObject('rgb I', [1000 1000 1000]));
light.addShape(pbrtShapeObject('disk', 'radius', 8));
light.addTransform(pbrtTransformObject('Translate', [ 0 9.9 0]));
light.addTransform(pbrtTransformObject('Rotate', [90 1 0 0]));
curPbrt.addLightSource(light);

frontOi = s3dRenderOIAndDepthMap(curPbrt, .050, []);
vcAddObject(frontOi); oiWindow;
%% flash rendering
% unix([fullfile(pbrtHome, '/src/bin/pbrt') ' ../cones/defaultBiggerZoom_Flash.pbrt']);
% % ISET will read the PBRT output
% % scene = sceneSet(scene,'fov', 8);
% oi = pbrt2oi('output_d.dat');
% % oi = oiSet (oi, 'horizontalfieldofview', 8 * 200/150 );
% vcAddAndSelectObject(oi);
% oiWindow;
% m = oiGet(oi, 'mean illuminance')

%% depth map rendering experiment
% This is now obsolete!

%basically, this renders a set number of depth maps, then combines them
%together.  They are combined using something that is NOT an averaging
%technique, such as using the median operation, or minimum operation.
%Initial evaluation shows that the median method works best.

%We plan to use rendertoolbox to make this rendering technique more
%streamlined and elegant in the future.
% 
% numDM =  31;  %number of depth map files to render
% 
% depthMap = zeros(300, 300, numDM);
% 
% for i = 1:numDM    
%     % depth map rendering
%     [path,name,ext] = fileparts(fname); 
%    
%     unix([fullfile(pbrtHome, '/src/bin/pbrt') path '/' name '_DM.pbrt' ]);
%     
%     dMapFile = 'output_d_dm_DM.dat'; 
%     depthMap(:,:,i) = s3dReadDepthMapFile(dMapFile);
% end
% 
% depthMapProcessedMin = min(depthMap, [], 3);
% depthMapProcessedMedian = median(depthMap, 3);
% figure; imagesc(depthMap(:,:,1));
% figure; imagesc(depthMapProcessedMin);
% figure; imagesc(depthMapProcessedMedian);
% 
% 
% %display the depth map
% %figure; imagesc(oi.depthMap);
% unix('cd ..');


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
