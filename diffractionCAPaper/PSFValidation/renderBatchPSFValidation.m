tic

%define template file
templateFile = 'pointTestTemplate';  %(no extension on here)
clear('templateFileArray');
templateFileArray{1} = batchFileClass(templateFile, '.pbrt');

%define conditions file
% conditionsFile = 'conditions1m.txt';
conditionsFile = 'conditionsPSFScratch.txt';
% conditionsFile = 'conditions2m50mmfvary.txt';

%change into correct directory
chdir(PSFValidationPath);
chdir('pointTest');

s3dRenderBatch(templateFileArray, conditionsFile);

toc
