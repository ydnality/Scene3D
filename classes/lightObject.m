% lightObject contains the superclass that makes 
classdef lightObject <  handle
    
    properties (SetAccess = private)
        name;
        type;
        spectrum; 
        propertyArray;
    end
    methods
        
        %default constructor.  The input variables may be omitted or left
        %with empty arguments if the user does not wish to specify them.  A
        %default value will be assumed.  
        function obj = lightObject(inName, inType, inSpectrum)
            if(ieNotDefined('inName'))
                obj.setName('defaultLight');
            else
                obj.setName(inName);
            end
            if(ieNotDefined('inType'))
                obj.setType('spot');
            else
                obj.setType(inType);
            end
            if(ieNotDefined('inSpectrum'))
                obj.setSpectrum(spectrumObject());
            else
                obj.setSpectrum(inSpectrum);
            end       
        end
        
        %sets the spectrum of the light.  inSpectrum must be of type
        %spectrumObject
        function setSpectrum(obj, inSpectrum)
           if (isa(inSpectrum, 'spectrumObject'))
                obj.spectrum = inSpectrum; 
           else
                error('**Warning: must be of type spectrumObject!!');
           end
           return;
        end

        function setType(obj, inType)
           obj.type = inType; 
           return;
        end

        function setName(obj, inName)
           obj.name = inName; 
           return;
        end    
        
        function writeFile(obj, fid)
            fprintf(fid,'\n\tLightSource\n');
            fprintf(fid,'\t\t"%s"\n', obj.type);
            fprintf(fid,'\t\t"%s" [', obj.spectrum.type);
            fprintf(fid,'%f ', obj.spectrum.value);
            fprintf(fid,']\n');
        end
    end
    
    methods (Abstract)
        move(obj, offset)
    end
end