%% Create dgauss50mm.mat

lensfileName = fullfile(s3dRootPath,'data','lens','dgauss.50mm.dat');
nSamples = 151;
apertureMiddleD = 10;   % mm
lens = lensMEObject('apertureSample', [nSamples nSamples], ...
    'fileName', lensFileName, ...
    'apertureMiddleD', apertureMiddleD);

lensfileName = fullfile(s3dRootPath,'data','lens','dgauss.50mm.mat');

save(lensfileName,'lens')

%%