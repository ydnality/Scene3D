%% loops through a series of different lights for a batch job
%makes a new tempPbrtFiles directory and puts all the files for the batch
%job in there.  Note that we must copy all the pbrt files in the current
%directory due to the way that s3dRenderOI is configured. 
filePath = [datapath '/validate/pbrtObject/']
chdir(filePath);
mkdir('batchPbrtFiles');
unix('rm batchPbrtFiles/*');
unix('cp * ./batchPbrtFiles/');
chdir('batchPbrtFiles');

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
   curPbrt.lightSourceArray{1}.setSpectrum(spectrumObject('rgb I', lightList(i, :))); %sets the spectrum to the one in teh list
   tmpFileName = ['deleteMe' int2str(i) '.pbrt']; %the pbrt file name
   curPbrt.writeFile(tmpFileName);
   oi = s3dRenderOI(tmpFileName, 50, [filePath '/batchPbrtFiles/'], tmpFileName);   %renders scene and displays as oi
end

chdir('..');
%the result should be the depthTargeSpheres scene with an additional green
%plane in the front.  The lighting should change between a series of
%different specifies lights.


%% this section changes the directions of the lights
%makes a new tempPbrtFiles directory and puts all the files for the batch
%job in there.  Note that we must copy all the pbrt files in the current
%directory due to the way that s3dRenderOI is configured. 
filePath = [datapath '/validate/pbrtObject/']
chdir(filePath);
mkdir('batchPbrtFiles');
unix('rm batchPbrtFiles/*');
unix('cp * ./batchPbrtFiles/');
chdir('batchPbrtFiles');

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

%loop through all light directions and render an image
for i = 1:size(lightOffsetList, 1)
   clear curPbrt;
   curPbrt = pbrtObject();
   curPbrt.lightSourceArray{1}.move(lightOffsetList(i, :));
   tmpFileName = ['deleteMe' int2str(i) '.pbrt']; %the pbrt file name
   curPbrt.writeFile(tmpFileName);
   oi = s3dRenderOI(tmpFileName, 50, [filePath '/batchPbrtFiles/'], tmpFileName);   %renders scene and displays as oi
end

chdir('..');
%the result should be the depthTargeSpheres scene with each one exhibiting
%slightly different lighting conditions.  The lights are moved around in a
%grid pattern.

