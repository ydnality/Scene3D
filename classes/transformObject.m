% wrapper class for a transform
classdef transformObject <  handle
%     enumeration
%        Translate ('Translate', [1 1 1])
%        Scale ('Scale', [1 1 1])
%     end
    properties (SetAccess = private)
        type;   %TODO: add enumeration for type
        data;
    end
    methods
        
        %default constructor.  The input variables may be omitted or left
        %with empty arguments if the user does not wish to specify them.  A
        %default value will be assumed.  
        function obj = transformObject(inType, inData)

            if(ieNotDefined('inType'))
                obj.setType('Translate');
            else
                obj.setType(inType);
            end
            if(ieNotDefined('inData'))
                obj.setData([0 0 0]);
            else
                obj.setData(inData);
            end       
        end

        %TODO: add error checking!
        function setType(obj, inType)
           validateattributes(inType, {'char'}, {'nonempty'});
           obj.type = inType; 
           return;
        end
        
        %TODO: error check for types
        function setData(obj, inData)
           validateattributes(inData, {'double'}, {'nonempty'});
           obj.data = inData; 
           return;
        end
        
        %write the pbrt output to file
        %TODO: check fid error checking
        function writeFile(obj, fid)
            fprintf(fid,'%s ', obj.type);
            fprintf(fid,'%f ', obj.data);
            fprintf(fid,'\n');
        end
    end
end