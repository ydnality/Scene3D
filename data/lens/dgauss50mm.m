%% Create dgauss.50mm.mat
%  11 element lens

lensFileDat = fullfile(s3dRootPath,'data','lens','dgauss.50mm.dat');
nSamples = 151;
apertureMiddleD = 10;   % mm
lens = lensMEObject('apertureSample', [nSamples nSamples], ...
    'fileName', lensFileDat, ...
    'apertureMiddleD', apertureMiddleD);

lens.draw

%%  Save
lensFileName = fullfile(s3dRootPath,'data','lens','dgauss.50mm.mat');
save(lensFileName,'lens')

%%  Load
clear lens
load(lensFileName,'lens')
lens.draw

%%