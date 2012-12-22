% Runs PBRT and imports it in ISET for a test point.  This one tries to
function output = s_s3dRenderChromaticAberrationTest(jobIndex)
% produce a high quality rendering, with chromatic abberation. 

%% PBRT will run the PBRT script

chdir(fullfile(s3dRootPath, 'scripts'));
outDir = ['tempDir' int2str(jobIndex)];
unix(['rm -rf ' outDir]);
unix(['mkdir ' outDir]);  %this forces the files to be generated here in scripts/testDir
chdir(outDir);

%to use this in proclus, we can simply make tons of new directories, and
%run the script in those directories.  use int2str to convert variables to string names  

% list of all chromatic aberration renderings - uncomment and run the one
% you wish to run

% slanted bar rendering
unix([fullfile(pbrtHome, '/src/bin/pbrt') ' ../pbrtFiles/chromaticAberration.pbrt']);

% radial lines rendering
%unix([fullfile(pbrtHome, '/src/bin/pbrt') ' chromaticAberrationRadial.pbrt']); 



%% ISET will read the PBRT output
% scene = sceneSet(scene,'fov', 8);
oi = pbrt2oi('chromaticAberration_d.dat');
% oi = oiSet (oi, 'horizontalfieldofview', 8 * 200/150 );
vcAddAndSelectObject(oi);
oiWindow;

m = oiGet(oi, 'mean illuminance')
unix('cd ..');
output = 0;  %indicates successful execution
end