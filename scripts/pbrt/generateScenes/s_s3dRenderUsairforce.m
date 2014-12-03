% Runs PBRT and imports it in ISET for the usairforce target.  This one tries to
% produce a high quality rendering, with diffraction. 

%% PBRT will run the PBRT script
chdir(fullfile(s3dRootPath, 'scripts', 'pbrtFiles'));
unix('/usr/share/pbrt-v2-spectral/src/bin/pbrt usairforceTest.pbrt'); 

%% ISET will read the PBRT output
% scene = sceneSet(scene,'fov', 8);
oi = pbrt2oi('usairforceTest_d.dat');
% oi = oiSet (oi, 'horizontalfieldofview', 8 * 200/150 );
vcAddAndSelectObject(oi);
oiWindow;

m = oiGet(oi, 'mean illuminance')
unix('cd ..');
