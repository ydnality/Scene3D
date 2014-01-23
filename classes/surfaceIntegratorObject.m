% surfaceIntegratorObject sets up the surface integrator
classdef surfaceIntegratorObject <  handle
    %     obj.surfaceIntegrator.type  = 'surfaceIntegrator';
    %     obj.surfaceIntegrator.surfIntType = 'directlighting';
    %     obj.surfaceIntegrator.maxdepth = 0;
    properties (SetAccess = private)
        type;
        maxDepth;
    end
    methods
        
        %default constructor
        function obj = surfaceIntegratorObject(inType, inMaxDepth)
            %sampler type
            if (ieNotDefined('inType'))
                obj.type = 'directlighting';
            else
                obj.type = inType;
            end
            
            if (ieNotDefined('inMaxDepth'))
                % Example lens
                obj.maxDepth = 0;
            else
                obj.maxDepth = inMaxDepth;
            end
        end
        
        %TODO: error checking
        function setType(obj, inType)
            obj.type = inType;
        end
        
        %TODO: error checking
        function setMaxDepth(obj, inMaxDepth)
            obj.maxDepth = inMaxDepth;
        end
        
        %TODO: might want to use a parameter list instead of hard coded max
        %depth
    end
end