%% loops through a series of different lights for a batch job
chdir([datapath '/tmp/']);

% spectrum list
lightList = ...
    [1000 0 0;
    0 1000 0;
    0 0 1000;
    1000 1000 1000];


%loop through all spectrums and render an image
for i = 1:size(lightList, 1)
   clear curPbrt;
   curPbrt = pbrtObject();
   curPbrt.lightSourceArray{1}.setSpectrum(spectrumObject('rgb I', lightList(i, :)));
   tmpFileName = ['deleteMe' int2str(i) '.pbrt']; %the pbrt file name
   curPbrt.writeFile(tmpFileName);
   oi = s3dRenderScene(tmpFileName, 50, [dataPath '/tmp/'], tmpFileName);   %renders scene and displays as oi
end

%the result should be the depthTargeSpheres scene with an additional green
%plane in the front.  The lighting should change between a series of
%different specifies lights.


%% this section changes the directions of the lights
chdir([datapath '/tmp/']);

% direction list
lightOffsetList = ...
    [0 0 0;
     1 0 0;
     2 0 0;
     0 0 1;
     1 0 1;
     2 0 1;
     0 0 2;
     1 0 2;
     2 0 2];

%loop through all spectrums and render an image
for i = 1:size(lightOffsetList, 1)
   clear curPbrt;
   curPbrt = pbrtObject();
   curPbrt.lightSourceArray{1}.move(lightOffsetList(i, :));
   tmpFileName = ['deleteMe' int2str(i) '.pbrt']; %the pbrt file name
   curPbrt.writeFile(tmpFileName);
   oi = s3dRenderScene(tmpFileName, 50, [dataPath '/tmp/'], tmpFileName);   %renders scene and displays as oi
end

%the result should be the depthTargeSpheres scene with each one exhibiting
%slightly different lighting conditions.  The lights are moved around in a
%grid pattern.

