% a simple pbrt geometry
classdef geometryObject <  handle
    
    %             examplePlane.type = 'geometry';
    %             examplePlane.name = 'Example Plane';
    %             examplePlane.material = 'greenLambertian';
    %             examplePlane.triangleMesh = [0 1 2;
    %                 0 2 3];
    %             examplePlane.points = [1.000000 1.000000 0.000000;
    %                 -1.000000 1.000000 0.000000;
    %                 -1.000000 -1.000000 0.000000;
    %                 1.000000 -1.000000 0.000000];
    %             examplePlane.transform = [-4.37113882867e-008 -0.0 -1.0 0.0;
    %                 1.0 -4.37113882867e-008 -4.37113882867e-008 0.0 ;
    %                 -4.37113882867e-008 -1.0 1.91068567692e-015 0.0;
    %                 0.0 0.0 8.87306690216 1.0];
    properties (SetAccess = private)
        name;
        material;
        triangleMesh;
        points;
        transform;
    end
    methods
        
        %default constructor
        function obj = geometryObject(inName, inMaterial, inTriMesh, inPoints, inTransform)
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
                obj.transform = [-4.37113882867e-008 -0.0 -1.0 0.0;
                    1.0 -4.37113882867e-008 -4.37113882867e-008 0.0 ;
                    -4.37113882867e-008 -1.0 1.91068567692e-015 0.0;
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
    end
    
end