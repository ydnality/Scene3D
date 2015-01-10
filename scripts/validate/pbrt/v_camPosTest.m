%% Loops through a series of different camera positions
%
% Requires PBRT on your path
%
% makes a new tempPbrtFiles directory and puts all the files for the batch
% job in there.  Note that we must copy all the pbrt files in the current
% directory due to the way that s3dRenderOI is configured.
%
% Andy to either delete or fix.
%
% Programming TODO
%
% Check whether pbrt exists gracefully
% Comment a little more
% Maybe we can use this for a unit test?  Save a simple scene and check
% that this thing always works?
%
%
%% camera position x-z offset list
cameraOffsetList = ...
    [0 0 0;
     1 0 0;
     2 0 0;
     0 0 1;
     1 0 1;
     2 0 1;
     0 0 2;
     1 0 2;
     2 0 2];

% geometry position list
geoList = [0 2 4 6];  %this is the x coordinate of the 

%loop through all camera positions in list and render
for i = 1:size(cameraOffsetList, 1)
   clear curPbrt;
   curPbrt = pbrtObject();
   
   %move camera
   curPbrt.camera.moveCamera(cameraOffsetList(i, :));   
   
   tmpFileName = ['deleteMe' int2str(i) '.pbrt'];
   curPbrt.writeFile(tmpFileName);
%    oi = s3dRenderOI(tmpFileName, 50, [filePath '/batchPbrtFiles/'], tmpFileName);
   oi = s3dRenderOI(tmpFileName, 50, tmpFileName);
end


%the result should be the depthTargeSpheres scene.  The camera will move
%left-right, and up-down depending on the offset specified.

%%