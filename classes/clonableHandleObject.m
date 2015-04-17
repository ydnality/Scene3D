classdef clonableHandleObject < handle
    %This is an abstract class (should not be used alone, but is meant to
    %be inherited only)
    %
    %The main functionality gained in addition to handle ability is the 
    %ability to clone an object (deep copy).
    %
    %code is inspired and/or copied from: 
    %http://www.mathworks.com/matlabcentral/fileexchange/22965-clone-handle-object-using-matlab-oop
    
    properties
    end
    
    methods
        
        function obj = makeDeepCopy(obj, oldObj)
            %clones the oldObj and makes a copy as the current object
            
             

            %only look at the current object for properties - this deals
            %with inheritance cases
            props = properties(obj);
            for i = 1:length(props)
                % Use Dynamic Expressions to copy the required property.
                % For more info on usage of Dylly namic Expressions, refer to
                % the section "Creating Field Names Dynamically" in:
                % web([docroot '/techdoc/matlab_prog/br04bw6-38.html#br1v5a9-1'])
%                 if isa(oldObj.(props{i}), 'clonableHandleObject')
%                     obj.(props{i}) = clonableHandleObject();
%                     obj.(props{i}).makeDeepCopy(oldObj.(props{i}));
%                 else
            	obj.(props{i}) = oldObj.(props{i});
                %if (iscell(props{i}))
                    
                %end
%                 end
            end
%             newObj = class(newObj, class(obj));   %TODO: FIX - want a
%             real copy that's automated
        end
        
    end
    
end