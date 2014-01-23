% a simple pbrt material
classdef materialObject <  handle
    
%     greenLambertian.type = 'material';
%     greenLambertian.name = 'greenLambertian';
%     greenLambertian.matType = 'matte';
%     greenLambertian.kd.type = 'kd';
%     greenLambertian.kd.kdType = 'color Kd';
%     greenLambertian.kd.value = [0 0.374624 0];
    
    properties (SetAccess = private)
        name;
        type;
        propertyArray;
    end
    methods
        
        %default constructor
        function obj = materialObject(inName, inType, inPropArray)
            if (ieNotDefined('inName'))
                obj.name = 'defaultMaterial';
            else
                obj.name = inName;
            end
            
            if (ieNotDefined('inType'))
                % Example lens
                obj.type = 'color Kd';
            else
                obj.type = inType;
            end
            
            if (ieNotDefined('inPropArray'))
                obj.propertyArray = cell(1,1);
                obj.propertyArray{1}.type = 'color Kd';  %TODO: make an actual object for this
                obj.propertyArray{1}.value = [0 0.374624 0];
            else
                obj.propertyArray = inPropArray;
            end
        end
        
        %TODO: error checking
        function setName(obj, inName)
            obj.name = inName;
        end
        
        %TODO: error checking
        function setType(obj, inType)
            obj.type = inType;
        end
        
        %TODO: error checking
        function addProperty(obj, inProperty)
            aLength = length(obj.propertyArray);
            obj.propertyArray{aLength+1} = inProperty;
        end
        
        %TODO: add a pop method? or some other array manipulation
    end
    
end