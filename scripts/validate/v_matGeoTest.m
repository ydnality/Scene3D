%% loops through a series of different surface and geometries (shifts the example plane in this example and changes its colors)
%makes a new tempPbrtFiles directory and puts all the files for the batch
%job in there.  Note that we must copy all the pbrt files in the current
%directory due to the way that s3dRenderOI is configured. 
% filePath = [datapath '/validate/pbrtObject/']
% chdir(filePath);
% mkdir('batchPbrtFiles');
% unix('rm batchPbrtFiles/*');
% unix('cp * ./batchPbrtFiles/');
% chdir('batchPbrtFiles');

% material list.  each line represents an rgb reflectance.
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
   newMaterial = pbrtMaterialObject(['mat' int2str(i)], 'matte', pbrtPropertyObject('color Kd', matList(i,:)));
   curPbrt.addMaterial(newMaterial); 
   
   %add new geoemtry
   %pbrtGeometryObject(inName, inMaterial, inTriMesh, inPoints, inTransform)
   %the plane is shifted in the x direction depending on the job number
   newTransform = [0 0 -1 0;
                    1 0 0 0 ;
                    0 -1 0 0;
                    geoList(i) 0 3  1]; %8.87306690216     %x direction is to the right, y is into the screen, z is up
   newGeometry = pbrtGeometryObject('newGeom', ['mat' int2str(i)], [], [], newTransform);  
   curPbrt.addGeometry(newGeometry);
   
   tmpFileName = ['deleteMe' int2str(i) '.pbrt'];
   curPbrt.writeFile(tmpFileName);
   oi = s3dRenderOI(curPbrt, 50, tmpFileName);
end

% chdir('..');

%the result should be the depthTargeSpheres scene with an additional plane
%on the bottom.  the plane should move to the right and change colors with
%subsequent jobs.