%% loops through a series of different surface and geometries

chdir([datapath '/tmp/']);

% camera position x-z offset list
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
   oi = s3dRenderScene(tmpFileName, 50, [dataPath '/tmp/'], tmpFileName);
end

%the result should be the depthTargeSpheres scene with an additional plane
%on the bottom.  the plane should move to the right and change colors with
%subsequent jobs.