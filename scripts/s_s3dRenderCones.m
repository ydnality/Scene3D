% Runs PBRT and imports it in ISET for the bench scene. 

%% PBRT will run the PBRT script
chdir(fullfile(s3dRootPath, 'scripts', 'pbrtFiles'));

% slanted bar rendering
unix([fullfile(pbrtHome, '/src/bin/pbrt') ' cones/defaultBiggerZoom.pbrt']);

%% ISET will read the PBRT output
% scene = sceneSet(scene,'fov', 8);
oi = pbrt2oi('cones_d.dat');
% oi = oiSet (oi, 'horizontalfieldofview', 8 * 200/150 );
vcAddAndSelectObject(oi);
oiWindow;

m = oiGet(oi, 'mean illuminance')
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
