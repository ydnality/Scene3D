% basic pbrt shapes such as spheres, cylinders etc.  Currently, only a
% sphere is supported
% TODO: use pbrtPropertyArrayObject instead
classdef pbrtShapeObject <  handle

    properties (SetAccess = private)
        type;   %TODO: add enumeration for type
        parameterArray;
        dataArray;
    end
    methods
        
        function obj = pbrtShapeObject(inType, inParameter, inData)
        %obj = pbrtShapeObject(inType, inParameter, inData)
        %
        %default constructor.  The input variables may be omitted or left
        %with empty arguments if the user does not wish to specify them.  A
        %default value will be assumed.  
        %inData and inParameter must both be specified, or none should be
        %specified.  
        %inParameter: one parameter that should be added to the parameter
        %array
        %inData: one data value corresponding to that parameter
        %TODO: support for multiple parameters or data
            if(ieNotDefined('inType'))
                obj.setType('sphere');
            else
                obj.setType(inType);
            end
            
            if(ieNotDefined('inParameter') || ieNotDefined('inData') )
                obj.addParameter('radius', .25);
            else
                obj.addParameter(inParameter, inData);
            end
        end

        function setType(obj, inType)
        %setType(obj, inType)
        %
        %sets the type of pbrtShape
        %
        %inType: must be of type character array (string)
           validateattributes(inType, {'char'}, {'nonempty'});
           obj.type = inType; 
           return;
        end
        
        function addParameter(obj, inParameter, inData)
        %addParameter(obj, inParameter, inData)
        %
        %inParameter: must be a character array (string).
        %inData: the input data.
        %TODO: error check for types
           validateattributes(inParameter, {'char'}, {'nonempty'});
           obj.parameterArray{length(obj.parameterArray) + 1} = inParameter; 
           
%            validateattributes(inData, {'double'}, {'nonempty'});
           validateattributes(inData, {'numeric'}, {'nonempty'});
           obj.dataArray{length(obj.dataArray) + 1} = inData; 
           return;
        end

        function writeFile(obj, fid)
        %writeFile(obj, fid)
        %
        %write the pbrt output to file
            fprintf(fid,'\tShape "%s" ', obj.type);
            
            for i = 1:length(obj.parameterArray)
                %TODO: add support for more shapes here...
                
                switch obj.parameterArray{i}
                    case 'radius'
                        fprintf(fid, '"float ');
                    case 'height'
                        fprintf(fid, '"float ');
                    otherwise
                        error(['unknown pbrtShapeObject parameter: ' obj.parameterArray{i}]);
                end
                %                 if (strcmp(obj.parameterArray{i}, 'radius'))
                %                     fprintf(fid, 'float ');
                %                 end
                fprintf(fid,'%s" [', obj.parameterArray{i});
                fprintf(fid,'%f ', obj.dataArray{i});
                fprintf(fid,']\n'); 
            end
        end
    end
end