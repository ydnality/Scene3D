%% PBRT will run the PBRT script
chdir(fullfile(s3dRootPath, 'scripts', 'validate'));
unix('/usr/share/pbrt-v2-spectral/src/bin/pbrt pointTest.pbrt'); 
%%IMPORTANT NOTE: please replace '/usr/share/pbrt-v2-spectral/src/bin/pbrt' with the location of your installation of pbrt

%% ISET will read the PBRT output
% scene = sceneSet(scene,'fov', 8);
oi = pbrt2oi('pointTest_d.dat');
% oi = oiSet (oi, 'horizontalfieldofview', 8 * 200/150 );
vcAddAndSelectObject(oi);
oiWindow;

m = oiGet(oi, 'mean illuminance')
