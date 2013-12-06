% Runs PBRT and imports it in ISET for a test point.  This one tries to
% produce a high quality rendering, with diffraction. 


chdir(PSFValidationPath);
fname = '../pointTest/pointTest25mm.pbrt'

mkdir('tempOutput');
chdir('tempOutput');
unix('rm *');

%% scene rendering
% unix([fullfile(pbrtHome, '/src/bin/pbrt') fname '--outfile output.dat']);
outfile = 'pointTest_out.dat';
unix([fullfile(pbrtHome, '/src/bin/pbrt ') fname ' --outfile ' outfile]);

% ISET will read the PBRT output
% scene = sceneSet(scene,'fov', 8);
oi = pbrt2oi(outfile);
% oi = oiSet (oi, 'horizontalfieldofview', 8 * 200/150 );
vcAddAndSelectObject(oi);
oiWindow;
m = oiGet(oi, 'mean illuminance')