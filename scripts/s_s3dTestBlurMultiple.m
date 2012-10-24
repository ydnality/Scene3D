%% load default table_sphere scene and depth map, depth map is pre-made
%% using separate function  

% Delete?

%TODO: figure out units!  - what does iset assume is size of optical image?
%TODO: aperture!  - how does this fit into thin lens equation?
%TODO: time of computation - right now very slow!


% note - w

% scene = s3dRT2Scene ([pwd '\image_data\picMat_orange_rad.mat'], 400:10:700); 
scene = s3dRT2Scene ([pwd '\image_data\picMat_orange_rad.mat'], 400:10:700); 
optics = opticsCreate('standard (1-inch)');
optics = opticsSet(optics , 'f#', 5.6);   % this does not affect anything right now - need to modify
optics = opticsSet(optics,'focallength',18e-3);  
scene = sceneSet(scene, 'samplesize', 2.6247e-003);
imgPlaneDist = 18.6e-3;

oiFinal = s3dRenderDepthDefocus(scene, optics, imgPlaneDist, 5);
% numIncrements = 5;
% s3dRenderDepthDefocus
vcAddAndSelectObject('oi',oiFinal); oiWindow;