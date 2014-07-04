% samplerObject contains the camera position, lens, sensor
classdef pbrtSamplerObject <  pbrtPropertyArrayObject
    %     obj.sampler.type = 'sampler';
    %     obj.sampler.samplerType = 'lowdiscrepancy';
    %     obj.sampler.pixelsamples = 512;
    properties (SetAccess = private)
        type;
%         pixelSamples;
    end
    methods
        
        %default constructor
        function obj = pbrtSamplerObject(inType, inProperty)
            %function obj = pbrtSamplerObject(inType, inProperty)
            
            %sampler type
            if (ieNotDefined('inType'))
                obj.type = 'lowdiscrepancy';
            else
                obj.type = inType;
            end
            
%             if (ieNotDefined('inPixelSamples'))
%                 % Example lens
%                 obj.pixelSamples = 128;
%             else
%                 obj.pixelSamples = inPixelSamples;
%             end
%           
            validateattributes(inProperty, {'pbrtPropertyObject'}, {'nonempty'});  
            obj.removeProperty();
            obj.addProperty(inProperty);
        end
        
        %TODO: error checking
        function setType(obj, inType)
            %function setType(obj, inType)
            obj.type = inType;
        end
        
%         %TODO: error checking
%         function setPixelSamples(obj, inPixelSamples)
%             obj.pixelSamples = inPixelSamples;
%         end
%         
        %TODO: might want to use a parameter list instead of hard coded
        %pixelSamples
        
        function writeFile(obj, fid)
            % function writeFile(obj, fid)
            
            fprintf(fid,'\n\nSampler "%s"\n', obj.type);    %TODO: consider putting this inside each of the objects? i'm not sure yet
            
            %assume numeric output for value for now
            %TODO: make this more general
            for i = 1:length(obj.propertyArray)
                if(ischar(obj.propertyArray{i}.value))
                    fprintf(fid,'\t"%s" %s\n',obj.propertyArray{i}.type, obj.propertyArray{i}.value);
                else
                    fprintf(fid,'\t"%s" [%i]\n',obj.propertyArray{i}.type, obj.propertyArray{i}.value);
                end
            end
%             fprintf(fid,'\t"integer pixelsamples" [%i]\n',obj.sampler.pixelSamples);
        end
    end
end