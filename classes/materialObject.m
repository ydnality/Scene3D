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
        %TODO: doccument
        function obj = materialObject(inName, inType, inProperty)
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
            
            obj.propertyArray = cell(1,1);
            if (ieNotDefined('inProperty'))
                obj.propertyArray{1} = propertyObject('color Kd', [0 0.374624 0]);
            else
                obj.propertyArray{1} = inProperty;
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
        
        %prints the pbrt representation to file object fid
        function returnVal = writeFile(obj, fid)
            fprintf(fid,'\n\nMakeNamedMaterial "%s"\n', obj.name);
            fprintf(fid,'\t"string type" ["%s"]\n', obj.type);
            
            %loop through property array and print the corresponding values
            for j = 1:length(obj.propertyArray)
                obj.propertyArray{j}.writeFile(fid);
            end
            returnVal = 1;
        end
    end
    
end