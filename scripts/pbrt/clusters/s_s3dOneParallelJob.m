% This function runs one single pbrt job on Proclus.  This is a wrapper function 
% for the pbrt call.  This funciton also adds a random pause between 0 and
% 1 seconds to every job, to allow for a proper random seed.
% inputs: pbrtFile: the name of the pbrt file to execute, 
% jobIndex: the current index for the job (used for random seed assignment)
% numpRrocs: the total number of processes used
function output = s_s3dOneParallelJob(pbrtFile, jobIndex, numProcs)
% produce a high quality rendering, with chromatic abberation. 
tic
%% PBRT will run the PBRT script
disp('jobIndex')
disp(int2str(jobIndex))
chdir(fullfile(s3dRootPath, 'scripts'));
outDir = ['tempDir' int2str(jobIndex)];
unix(['rm -rf ' outDir]);
unix(['mkdir ' outDir]);  %this forces the files to be generated here in scripts/testDir
chdir(outDir);

%a problem with PBRT is that the random seed is assigned using the time
%if these processes are submitted simultaneously, all random starting points
%will be the same - so all images are identical - this is NOT desired
%Thus, we make a random delay period over here
randomState = rng(jobIndex); %sets the seed based off of job index
rng(randomState); %initializes random number based off of seed
%for some unknown reason - the above attempt does not work, but here is a workaround
randomPause = rand(1,numProcs);
randomPause = randomPause(jobIndex);
disp(sprintf('pausing for %f seconds.', randomPause)); 
pause(randomPause);


%to use this in proclus, we can simply make tons of new directories, and
%run the script in those directories.  use int2str to convert variables to string names  

% list of all chromatic aberration renderings - uncomment and run the one
% you wish to run

% slanted bar rendering
unix([fullfile(pbrtHome, '/src/bin/pbrt ..') pbrtFile]);

% radial lines rendering
%unix([fullfile(pbrtHome, '/src/bin/pbrt') ' chromaticAberrationRadial.pbrt']); 



%% ISET will read the PBRT output
% scene = sceneSet(scene,'fov', 8);
% oi = pbrt2oi('chromaticAberration_d.dat');
% oi = oiSet (oi, 'horizontalfieldofview', 8 * 200/150 );
%vcAddAndSelectObject(oi);
%oiWindow;

m = oiGet(oi, 'mean illuminance')
unix('cd ..');
output = 0;  %indicates successful execution
toc
end
