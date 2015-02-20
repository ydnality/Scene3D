% pbrtLightObject contains the superclass that makes 
classdef pbrtLightObject <  pbrtPropertyArrayObject
    
    properties (SetAccess = private)
        name;
        type;
        spectrum; 
        %propertyArray;  
    end
    methods
        
        %default constructor.  The input variables may be omitted or left
        %with empty arguments if the user does not wish to specify them.  A
        %default value will be assumed.  
        function obj = pbrtLightObject(inName, inType, inSpectrum)
            if(ieNotDefined('inName'))
                obj.setName('defaultLight');
            else
                obj.setName(inName);
            end
            if(ieNotDefined('inType'))
                obj.setType('point');
            else
                obj.setType(inType);
            end
            if(ieNotDefined('inSpectrum'))
                obj.setSpectrum(pbrtSpectrumObject());
            else
                obj.setSpectrum(inSpectrum);
            end       
%             if(~ieNotDefined('inPropertyArray'))
%                 obj.propertyArray = inPropertyArray;
%             end
        end
        
        %sets the spectrum of the light.  inSpectrum must be of type
        %pbrtSpectrumObject
        function setSpectrum(obj, inSpectrum)
           if (isa(inSpectrum, 'pbrtSpectrumObject'))
                obj.spectrum = inSpectrum; 
           else
                error('**Warning: must be of type pbrtSpectrumObject!!');
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
    
    %methods (Abstract)
    %    move(obj, offset)
    %end
end