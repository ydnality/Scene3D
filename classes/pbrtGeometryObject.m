% a simple pbrt geometry
classdef pbrtGeometryObject <  handle
    properties (SetAccess = private)
        name;
        material;
        triangleMesh;
        points;
        transform;
    end
    methods
        
        %default constructor
        function obj = pbrtGeometryObject(inName, inMaterial, inTriMesh, inPoints, inTransform)
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
        
        %TODO: error checking
        function setName(obj, inName)
            obj.name = inName;
        end
        
        %TODO: error checking
        function setMaterial(obj, inMat)
            obj.material = inMat;
        end
        
        %TODO: error checking
        function setTransform(obj, inTransform)
            obj.transform = inTransform;
        end
        
        %TODO: error checking
        function setPoints(obj, inPoints)
            obj.points = inPoints;
        end

        %TODO: error checking
        function setTriangleMesh(obj, inTriangleMesh)
            obj.triangleMesh = inTriangleMesh;
        end
        
        %prints the pbrt representation to file object fid
        function returnVal = writeFile(obj, fid)
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