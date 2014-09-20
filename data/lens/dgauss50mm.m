%% Create dgauss.50mm.mat
%  11 element lens

lensFileDat = fullfile(s3dRootPath,'data','lens','dgauss.50mm.dat');
nSamples = 151;
apertureMiddleD = 10;   % mm
lens = lensC('apertureSample', [nSamples nSamples], ...
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

%% Now create dgauss.50mm.2.mat

% This is the same except the final aperture is 2mm.  AL finds this useful
% for some PBRT calculations.
lensFileDat = fullfile(s3dRootPath,'data','lens','dgauss.50mm.2mm.dat');
nSamples = 151;
apertureMiddleD = 10;   % mm
lens = lensC('apertureSample', [nSamples nSamples], ...
    'fileName', lensFileDat, ...
    'apertureMiddleD', apertureMiddleD);

lens.draw

%%  Save
lensFileName = fullfile(s3dRootPath,'data','lens','dgauss.50mm.2mm.mat');
save(lensFileName,'lens')

%%  Load
clear lens
load(lensFileName,'lens')
lens.draw
