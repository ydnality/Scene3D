% propertyObject
% TODO: document what this does
classdef propertyObject <  handle
    properties 
        type;
        value;
    end
    methods
        
        %default constructor
        function obj = propertyObject(inType, inValue)
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

    end
end