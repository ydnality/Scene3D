%% loops through a series of different lights for a batch job

clear templatePbrt;
templatePbrt = pbrtObject();

chdir([datapath '/tmp/']);

% spectrum list
lightList = ...
    [1000 0 0;
    0 1000 0;
    0 0 1000;
    1000 1000 1000];


%loop through all spectrums and render an image
for i = 1:size(lightList, 1)
   curPbrt = templatePbrt;
   curPbrt.lightSourceArray{1}.setSpectrum(spectrumObject('rgb I', lightList(i, :)));
   tmpFileName = ['deleteMe' int2str(i) '.pbrt'];
   curPbrt.writeFile(tmpFileName);
   oi = s3dRenderScene(tmpFileName, 50, [dataPath '/tmp/'], tmpFileName);
end

