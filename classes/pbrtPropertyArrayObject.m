% a class that is meant to be inherited that handles the propertyArray
classdef pbrtPropertyArrayObject <  handle
    properties (SetAccess = private)
        propertyArray;
    end
    methods
        
        %default constructor
        %TODO: doccument
        function obj = pbrtPropertyArrayObject(inProperty)
            obj.propertyArray = cell(1,1);
            if (ieNotDefined('inProperty'))  %default property
                obj.propertyArray{1} = pbrtPropertyObject('color Kd', [0 0.374624 0]);
            else
                obj.propertyArray{1} = inProperty;
            end
        end
        
        
        %TODO: error checking
        function addProperty(obj, inProperty)
            validateattributes(inProperty, {'pbrtPropertyObject'}, {'nonempty'});
            obj.propertyArray{end+1} = inProperty;
        end
 

        %removes the shape corresponding to the specified index
        %if deleteIndex is undefined, remove from the end
        %returns the deleted value
        function returnVal = removeProperty(obj, removeIndex)
            if (ieNotDefined('removeIndex'))
                returnVal = obj.propertyArray{end};
                obj.propertyArray(end)= [];
            else
                 validateattributes(removeIndex, {'numeric'}, {'positive'});
                returnVal = obj.propertyArray{removeIndex};
                obj.propertyArray(removeIndex)= [];
            end
        end
        
        %TODO: add a pop method? or some other array manipulation
    end
    
end