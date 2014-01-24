%% loops through a series of different surface and geometries

chdir([datapath '/tmp/']);

% material list
matList = ...
    [1 1 1;
    1 0 0;
    0 1 0;
    0 0 1];

% geometry position list
 geoList = [0 2 4 6];  %this is the x coordinate of the 

%loop through all materials and geometries and render an image
for i = 1:size(matList, 1)
   clear curPbrt;
   curPbrt = pbrtObject();
   %add a new material
   newMaterial = materialObject(['mat' int2str(i)], 'matte', propertyObject('color Kd', matList(i,:)));
   curPbrt.addMaterial(newMaterial); 
   
   %add new geoemtry
   %geometryObject(inName, inMaterial, inTriMesh, inPoints, inTransform)
   newTransform = [0 0 -1 0;
                    1 0 0 0 ;
                    0 -1 0 0;
                    geoList(i) 0 3  1]; %8.87306690216     %x direction is to the right, y is into the screen, z is up
   newGeometry = geometryObject('newGeom', ['mat' int2str(i)], [], [], newTransform);  
   curPbrt.addGeometry(newGeometry);
   
   tmpFileName = ['deleteMe' int2str(i) '.pbrt'];
   curPbrt.writeFile(tmpFileName);
   oi = s3dRenderScene(tmpFileName, 50, [dataPath '/tmp/'], tmpFileName);
end

%the result should be the depthTargeSpheres scene with an additional plane
%on the bottom.  the plane should move to the right and change colors with
%subsequent jobs.