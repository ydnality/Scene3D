%% loops through a series of different shape geometries (shifts the example white sphere in this example)
%makes a new tempPbrtFiles directory and puts all the files for the batch
%job in there. 

% geometry position list
geoList = [0 2 4 6];  %this is the x coordinate of the 

%loop through all materials and geometries and render an image
for i = 1:length(geoList)
   clear curPbrt;
   curPbrt = pbrtObject();
   %add a new material
   newMaterial = pbrtMaterialObject('mat1', 'matte', pbrtPropertyObject('color Kd', [1 1 1]));
   curPbrt.addMaterial(newMaterial); 
   
   %add new geoemtry
   %pbrtGeometryObject(inName, inMaterial, inTriMesh, inPoints, inTransform)
   %the plane is shifted in the x direction depending on the job number
   newTransform = [0 0 -1 0;
                    1 0 0 0 ;
                    0 -1 0 0;
                    geoList(i) 0 3  1]; %8.87306690216     %x direction is to the right, y is into the screen, z is up
   newGeometry = pbrtGeometryObject('sphere', 'mat1', pbrtShapeObject('sphere', 'radius', 2), [], newTransform);  
   curPbrt.addGeometry(newGeometry);
   
   tmpFileName = ['deleteMe' int2str(i) '.pbrt'];
   curPbrt.writeFile(tmpFileName);
   oi = s3dRenderOI(curPbrt, 50, tmpFileName);
end

%the result should be the depthTargeSpheres scene with an additional sphere
%on the bottom.  the plane should move to the right.