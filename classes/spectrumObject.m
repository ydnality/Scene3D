% lightObject contains the object to create a light in PBRT
classdef spectrumObject <  handle 
    
    properties %(SetAccess = private)
        type;
        value;
    end
    methods
        
        %default constructor
        function obj = spectrumObject(inType, inValue)
            if (ieNotDefined('inType'))
                obj.type = 'rgb I';
            else
                obj.type = inType;
            end
            
            if (ieNotDefined('inValue'))
                obj.value = [1000 1000 1000];
            else
                obj.value = inValue;
            end
        end
        
        function setType(obj, inType)
            obj.type = inType;
        end
        
        function setValue(obj, inValue)
            obj.value = inValue;
        end
    end 
end