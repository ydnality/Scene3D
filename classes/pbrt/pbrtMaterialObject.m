% A pbrt material.
classdef pbrtMaterialObject <  handle
    
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
       
        function obj = pbrtMaterialObject(inName, inType, inProperty)
        %obj = pbrtMaterialObject(inName, inType, inProperty)
        %inName: name of the material (default: 'defaultMaterial')
        %inType: type of material (default: 'color Kd')
        %inProperty: material property (default: pbrtPropertyObject('color
        %Kd', [0 0.374624 0]))
        
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
                obj.propertyArray{1} = pbrtPropertyObject('color Kd', [0 0.374624 0]);
            else
                obj.propertyArray{1} = inProperty;
            end
        end
        
        function setName(obj, inName)
        %setName(obj, inName)
        %sets the name of the material
        %
        %inName: Name to set the material as.  Must be a string.
        %TODO: error checking
            obj.name = inName;
        end
        
        function setType(obj, inType)
        %setType(obj, inType)
        %sets the type of the material
        %
        %inName: Type to set the material as.  Must be a string.
        %TODO: error checking
            obj.type = inType;
        end
        
        function addProperty(obj, inProperty)
        %addProperty(obj, inProperty)
        %adds a material property
        %
        %inProperty: the added property.  Must be a pbrtPropertyObject.   
        %TODO: error checking
        %TODO: add a pop method? or some other array manipulation
            aLength = length(obj.propertyArray);
            obj.propertyArray{aLength+1} = inProperty;
        end
        
        function returnVal = writeFile(obj, fid)
        %returnVal = writeFile(obj, fid)
        %prints the pbrt representation to file object fid
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