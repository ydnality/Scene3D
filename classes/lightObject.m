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
        
        %default constructor
        function obj = lightObject()
            obj.name = 'defaultLight';
            obj.type = 'spot';
            obj.spectrum = spectrumObject();
            obj.coneAngle = 180;
            obj.coneDeltaAngle = 180;
            obj.from =  [4.5 -90 8.5];
            obj.to = [4.5 -89 8.5];
        end
        
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
           obj.from = inFrom;
           obj.to = inTo;
           return;
        end
        
         %sets the direction of the camera
        %inFrom: how much the light should be moved. Note that the
        %direction stays the same
        function move(obj, offset)
           obj.from = obj.from + offset;
           obj.to = obj.to + offset;
           return;
        end
               
    end
end