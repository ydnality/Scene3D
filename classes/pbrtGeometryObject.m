% A simple pbrt geometry
classdef pbrtGeometryObject <  handle
    properties (SetAccess = private)
        name;
        material;
        triangleMesh;
        points;
        transform;
    end
    methods

        function obj = pbrtGeometryObject(inName, inMaterial, inTriMesh, inPoints, inTransform)
        %obj = pbrtGeometryObject(inName, inMaterial, inTriMesh, inPoints, inTransform)
        %default
        %
        %inName: name of the geometry (default: 'defaultGeometry').  Must
        %be a string.
        %inMaterial: name of the material to use for this geometry
        %(default: 'greenLambertian').  Must be a string.
        %inTriMesh: the triangle mesh that comprises the geometry (default: 
        %[0 1 2; 0 2 3];).  
        %inPoints: the points in space that correspond to the triangle mesh
        %alias numbers (default: [1.000000 1.000000 0.000000;
        %            -1.000000 1.000000 0.000000;
        %            -1.000000 -1.000000 0.000000;
        %            1.000000 -1.000000 0.000000];).  
        %inTransform: the transform of the geometry (4x4).  
        %           i.e.    [0 0 -1 0;
        %                   1 0 0 0 ;
        %                   0 -1 0 0;
        %                   0 0 0  1];
        %TODO: supply more examples.  
        
            if (ieNotDefined('inName'))
                obj.name = 'defaultGeometry';
            else
                obj.name = inName;
            end
            
            if (ieNotDefined('inMaterial'))
                obj.material = 'greenLambertian';  %TODO: change default later
            else
                obj.material = inMaterial;
            end
            
            %TODO:error checking
            if (ieNotDefined('inTriMesh'))
                obj.triangleMesh = [0 1 2; 0 2 3];
            else
                obj.triangleMesh = inTriMesh;
            end
            
            %TODO:error checking
            if (ieNotDefined('inPoints'))
                obj.points = [1.000000 1.000000 0.000000;
                    -1.000000 1.000000 0.000000;
                    -1.000000 -1.000000 0.000000;
                    1.000000 -1.000000 0.000000];
            else
                obj.points = inPoints;
            end
            
            %TODO:error checking
            if (ieNotDefined('inTransform'))
                obj.transform = [0 -0.0 -1.0 0.0;
                    1.0 -0 0 0.0 ;
                    0 -1.0 0 0.0;
                    0.0 0.0 8.87306690216 1.0];
            else
                obj.transform = inTransform;
            end            
        end

        function setName(obj, inName)
        %setName(obj, inName)
        %
        %Sets the name.
        %inName: must be a string.
        %TODO: error checking
            obj.name = inName;
        end
        
        
        function setMaterial(obj, inMat)
        %setMaterial(obj, inMat)
        %
        %Sets the material.
        %inMat: must be a string.  
            %TODO: error checking
            obj.material = inMat;
        end
        
        function setTransform(obj, inTransform)
        %setTransform(obj, inTransform)
        %
        %Sets the transform.
        %inPoints: must be 4x4 matrix
        %TODO: error checking
            obj.transform = inTransform;
        end
        
        function setPoints(obj, inPoints)
        %setPoints(obj, inPoints)
        %
        %Sets the points.
        %inPoints: must be a vector of floats of size n x 3, where n is the
        %number of points.
        %TODO: error checking
            obj.points = inPoints;
        end
        
        function setTriangleMesh(obj, inTriangleMesh)
        %setTriangleMesh(obj, inTriangleMesh)
        %
        %Sets the triangle mesh.
        %inTriangleMesh: must contain a matrix of integers of size n x 3,
        %where n is the number of triangles. 
        %TODO: error checking
            obj.triangleMesh = inTriangleMesh;
        end
        
        function returnVal = writeFile(obj, fid)
        %prints the pbrt representation to file object fid
            fprintf(fid,'\n\nAttributeBegin #%s\n', obj.name);
            
            fprintf(fid,'\n\tTransform \n\t[\n');
            fprintf(fid,'\t%f %f %f %f \n', obj.transform' );
            fprintf(fid,'\t]\n');
            
            fprintf(fid,'\tNamedMaterial "%s"\n', obj.material);
            
            fprintf(fid,'\tShape "trianglemesh" "integer indices" \n\t[\n');
            fprintf(fid,'\t%i %i %i\n', obj.triangleMesh');
            fprintf(fid,'\t]\n');
            
            fprintf(fid,'\t"point P" \n\t[\n');
            fprintf(fid,'\t%f %f %f\n', obj.points');
            fprintf(fid,'\t]\n');
            
            fprintf(fid,'\n\nAttributeEnd\n');
            returnVal = 1;
        end  
    end
end