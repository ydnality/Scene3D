
%for some reason - does NOT work with desktop mode
%% this section is to set up and run the parallel jobs
n = 10;
numProc = 24;
sgerun2([],['testparallel' int2str(n)], 0, 1:numProc);
s_s3dRenderChromaticAberrationTest(jobindex, 24) %need a more elegant solution later... constants maybe?
.

%% this line is to view the Matlab output
unix(['cat ~/sgeoutput/job_testparallel' int2str(n) '.o*'])
