% propertyObject - contains the type and value of each property.  This is
% meant to be put in array for general properties.
% TODO: document what this does
classdef pbrtPropertyObject <  handle
    properties 
        type;
        value;
    end
    methods
        
        %default constructor
        function obj = pbrtPropertyObject(inType, inValue)
            %sampler type
            if (ieNotDefined('inType'))
                obj.type = 'color Kd';
            else
                obj.type = inType;
            end
            
            if (ieNotDefined('inValue'))
                % Example lens
                obj.value = [1 1 1];
            else
                obj.value = inValue;
            end
        end
        
        %prints the pbrt representation to file object fid
        function returnVal = writeFile(obj, fid)
            fprintf(fid,'\t"%s" [', obj.type);
            fprintf(fid,'%f ', obj.value);
            fprintf(fid,']\n');
            returnVal = 1;
        end

    end
end