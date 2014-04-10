classdef pbrtSpectrumObject < handle
%defines a spectrum for PBRT.  
%the default input is of type 'rgb I' and value [1000 1000 1000].
%See PBRT documentation for additional inputs. 
%TODO: add more input types for users
    properties %(SetAccess = private)
        type;
        value;
    end
    methods
        
        %default constructor
        function obj = pbrtSpectrumObject(inType, inValue)
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