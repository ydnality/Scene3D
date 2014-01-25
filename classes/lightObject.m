% lightObject contains the object to create a light in PBRT
classdef lightObject <  handle
    
    properties (SetAccess = private)
        name;
        type;
        spectrum; 
        coneAngle;
        coneDeltaAngle;
        from;
        to;
    end
    methods
        
        %default constructor.  The input variables may be omitted or left
        %with empty arguments if the user does not wish to specify them.  A
        %default value will be assumed.  
        function obj = lightObject(inName, inType, inSpectrum, inConeAngle, inDeltaAngle, inFrom, inTo)
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
            if(ieNotDefined('inConeAngle'))
                obj.setAngle(180);
            else
                obj.setAngle(inConeAngle);
            end
            if(ieNotDefined('inDeltaAngle'))
                obj.setDeltaAngle(180);
            else
                obj.setDeltaAngle(inDeltaAngle);
            end
            if(ieNotDefined('inFrom'))
                obj.setFrom([4.5 -90 8.5]);
            else
                obj.setFrom(inFrom);
            end
            if(ieNotDefined('inTo'))
                obj.setTo([4.5 -89 8.5]);
            else
                obj.setTo(inTo);
            end            
        end
        
        %sets the spectrum of the light.  inSpectrum must be of type
        %spectrumObject
        function setSpectrum(obj, inSpectrum)
           if (isa(inSpectrum, 'spectrumObject'))
                obj.spectrum = inSpectrum; 
           else
                disp('**Warning: must be of type spectrumObject!!');
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
                
        function setAngle(obj, inConeAngle)
           obj.coneAngle = inConeAngle; 
           return;
        end    
        
        function setDeltaAngle(obj, inDeltAngle)
           obj.coneDeltaAngle = inDeltAngle; 
           return;
        end    
        
        %sets the direction of the camera
        %inFrom: will be the input initial point of vector representing
        %direction
        %inTo: is the second point of the vector
        function setDirection(obj, inFrom, inTo)
           obj.setFrom(inFrom);
           obj.setTo(inTo);
           return;
        end

        %sets the direction of the camera
        %inTo: is the second point of the vector of the direction
        function setTo(obj, inTo)
           obj.to = inTo;
           return;
        end        
        
        %sets the direction of the camera
        %inTo: is the first point of the vector of the direction
        function setFrom(obj, inFrom)
           obj.from = inFrom;
           return;
        end        
        
        %sets the direction of the camera
        %inFrom: how much the light should be moved. Note that the
        %direction stays the same
        function move(obj, offset)
           obj.setFrom(obj.from + offset);
           obj.setTo(obj.to + offset);
           return;
        end
               
    end
end