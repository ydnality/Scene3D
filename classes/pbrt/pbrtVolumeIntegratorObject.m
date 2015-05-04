% surfaceIntegratorObject sets up the surface integrator
classdef pbrtVolumeIntegratorObject <  handle
    %     obj.surfaceIntegrator.type  = 'surfaceIntegrator';
    %     obj.surfaceIntegrator.surfIntType = 'directlighting';
    %     obj.surfaceIntegrator.maxdepth = 0;
    properties (SetAccess = private)
        type;
        stepSize;
    end
    methods
        
        %default constructor
        function obj = pbrtVolumeIntegratorObject(inType, inStepSize)
            %sampler type
            if (ieNotDefined('inType'))
                obj.type = 'emission';
            else
                obj.type = inType;
            end
            
            if (ieNotDefined('inStepSize'))
                % Example lens
                obj.stepSize = 16;
            else
                obj.stepSize = inStepSize;
            end
        end
        
        %TODO: error checking
        function setType(obj, inType)
            obj.type = inType;
        end
        
        %TODO: error checking
        function setStepSize(obj, inStepSize)
            obj.stepSize = inStepSize;
        end
        
        %TODO: might want to use a parameter list instead of hard coded max
        %depth
        
        %prints the pbrt representation to file object fid
        function returnVal = writeFile(obj, fid)
            fprintf(fid,'\n\nVolumeIntegrator "%s"\n', obj.type);
            fprintf(fid,'\t"float stepsize" [%i]\n',obj.stepSize);
            returnVal = 1;
        end
        
    end
end