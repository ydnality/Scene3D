% cameraObject contains the camera position, lens, sensor
classdef samplerObject <  handle
    %     obj.sampler.type = 'sampler';
    %     obj.sampler.samplerType = 'lowdiscrepancy';
    %     obj.sampler.pixelsamples = 512;
    properties (SetAccess = private)
        type;
        pixelSamples;
    end
    methods
        
        %default constructor
        function obj = samplerObject(inType, inPixelSamples)
            %sampler type
            if (ieNotDefined('inType'))
                obj.type = 'lowdiscrepancy';
            else
                obj.type = inType;
            end
            
            if (ieNotDefined('inPixelSamples'))
                % Example lens
                obj.pixelSamples = 128;
            else
                obj.pixelSamples = inPixelSamples;
            end
        end
        
        %TODO: error checking
        function setType(obj, inType)
            obj.type = inType;
        end
        
        %TODO: error checking
        function setPixelSamples(obj, inPixelSamples)
            obj.pixelSamples = inPixelSamples;
        end
        
        %TODO: might want to use a parameter list instead of hard coded
        %pixelSamples
    end
end