%This script combines sub images created by the cluster.  
%numProc must be modified to reflect the correct number of subprocesses
%used.
%Note: must run this from the scripts directory

tic
numProc = 128;  %keep track of the number of parallel processes, and modify this number here

chdir('tempDir1');  %don't change this
outputFilename = 'output_d.dat';
%outputFilename = 'cones_d.dat';

firstOi = pbrt2oi(outputFilename);  %please keep track of the output file name, and place it here
vcAddAndSelectObject(firstOi);
cd ..
%for some reason pbrt2oi does a cd.. already - should change this...
tempPhotons = oiGet(firstOi, 'photons');

for i = 2:numProc
          chdir(['tempDir' int2str(i)])
oi = pbrt2oi(outputFilename);
tempPhotons = tempPhotons + oiGet(oi, 'photons');
    %vcAddAndSelectObject(oi); 
    cd ..
end

tempPhotons = tempPhotons./numProc;
newOi = oi;
newOi = oiSet(newOi, 'photons', tempPhotons);
vcAddAndSelectObject(newOi); oiWindow;
toc
'approximate time needed to combine images'