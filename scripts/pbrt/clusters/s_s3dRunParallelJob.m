%This script sets up and runs a job on the Proclus cluster
%jobName, and numProc, and pbrtFile can all be modified to reflect the
%desired configuration
%Note:sgerun2 is a cluster-specific matlab function written by Kendrick Kay,
%meant to load multiple parallel jobs on the cluster

%% initialize - only run this once
n = 1;

%% this section is to set up and run the parallel jobs - you can keep running this section for debugging
tic
%you may change these variables to fit your job
numProc = 128;  %number of parallel processes
jobName = 'bench';
n = n + 1; %this is the job ID - increment it when adding a new job
pbrtFile = '/pbrtFiles/chromaticAberration.pbrt'
%pbrtFile = '/pbrtFiles/benchScene/defaultBiggerHQP.pbrt'
%pbrtFile = '/pbrtFiles/benchScene/defaultBiggerCameraArray.pbrt'
%pbrtFile = '/pbrtFiles/cones/defaultBiggerZoom.pbrt'
%pbrtFile = '/pbrtFiles/bunnyScene/bunnies.pbrt'

sgerun2(['s_s3dOneParallelJob(''' pbrtFile ''', jobindex,'  int2str(numProc) ');'],[jobName int2str(n)], 0, 1:numProc);
toc
'approximate time needed to set up jobs'

'current time'
c = clock;
'hour'
c(4)
'minute' 
c(5)
%% this line is to view the Matlab output
unix(['cat ~/sgeoutput/job_' jobName int2str(n) '.o*'])
