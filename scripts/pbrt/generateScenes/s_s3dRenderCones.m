% Runs PBRT and imports it in ISET for the cones scene. 
%
% Creates the light field OI that we use to illustrate transverse chromatic
% aberration.  
%
% The number of samples needs to be increased to make a reasonable image.
% This could be made more efficient using MP's code based on the effective
% aperture.  In general, if we could integrate the bbm code with what we
% are doing, that would be good.
%
% AL Copyright Vistasoft Team 2015

%%
ieInit;

%% This section runs the PBRT script directly (demo) 
% dockerFlag = false;  % Need to have a local version of PBRT
% fname = fullfile(dataPath, 'pbrtScenes', 'cones', 'defaultBiggerZoom.pbrt');  %defaultBiggerZoom_depth.pbrt is for depth map
% oi = s3dRenderOIAndDepthMap(fname, 'transAb', dockerFlag);
% vcAddObject(oi); oiWindow;

%% This section uses a PBRT wrapper for more flexibility for rendering

%initialization
clear curPbrt;
curPbrt = pbrtObject();

includeFile = fullfile(dataPath, 'pbrtScenes', 'cones', 'cones.pbrt');
if ~exist(includeFile,'file'), error('Missing key include file'); end

curPbrt.addInclude(includeFile);
curPbrt.camera.addTransform(pbrtTransformObject('Scale', [1000 1000 1000]));
curPbrt.camera.addTransform(pbrtTransformObject('Rotate',  [-3 1 0 0]));
curPbrt.camera.addTransform(pbrtTransformObject('Rotate', [52 0 1 0]));
curPbrt.camera.addTransform(pbrtTransformObject('Translate', [-2.3 -.05 .5]));

%% Make a lens and add to pbrt object

filmDist = 40;  % This brings the front part of the image into reasonable focus
filmDiag = 10;
nSamples = 2048;  % This is too few, but doesn't take too long to run.
specFile = 'dgauss.50mm.dat';
apertureDiameter    = 3;  % mm
diffraction         = false;
chromaticAberration = true;

% set sampler
curPbrt.sampler.removeProperty();
curPbrt.sampler.addProperty(pbrtPropertyObject('integer pixelsamples', nSamples));

lens = pbrtLensRealisticObject(filmDist, filmDiag, specFile, apertureDiameter, diffraction, chromaticAberration, []);
curPbrt.camera.setLens(lens);

%% make an area light and add to pbrt object
light = pbrtAreaLightObject('area', pbrtSpectrumObject('color L', [1000 1000 1000]));
areaLightShape = pbrtShapeObject('disk', 'radius', 8);
light.addShape(areaLightShape);
light.addTransform(pbrtTransformObject('Translate', [ 0 9.9 0]));
light.addTransform(pbrtTransformObject('Rotate', [90 1 0 0]));

% light.removeProperty();
light.addProperty(pbrtPropertyObject('integer nsamples', 4));

% Remove the default and add the one that was just made
curPbrt.removeLight();
curPbrt.addLightSource(light);

% An optional alternative light source
% curPbrt.removeLight();
% curPbrt.addLightSource(pbrtLightInfiniteObject('infiniteLight', 16, [], [], []));
% 

% Remove the defaults
curPbrt.removeMaterial();
curPbrt.removeGeometry();

%% Run it

oiName = 'coneArray';
dockerFlag = true;
frontOi = s3dRenderOIAndDepthMap(curPbrt, oiName, dockerFlag);

% Visualize it
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
