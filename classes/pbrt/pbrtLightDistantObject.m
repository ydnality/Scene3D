% spotlightObject contains the subclass that makes pbrt spotlights
classdef pbrtLightDistantObject <  pbrtLightObject
    
    properties (SetAccess = private)
        from;
        to;
    end
    methods
        
        function obj = pbrtLightDistantObject(inName, inSpectrum, inFrom, inTo)
            %obj = pbrtLightSpotObject(inName, inSpectrum, inConeAngle, inDeltaAngle, inFrom, inTo)
            %default constructor.  The input variables may be omitted or left
            %with empty arguments if the user does not wish to specify them.  A
            %default value will be assumed.
        
            %superclass properties
            obj@pbrtLightObject();
            obj.setType('distant');
            
            if (~ieNotDefined('inName'))
                obj.setName(inName);
            end
            if (~ieNotDefined('inSpectrum'))
                obj.setSpectrum(inSpectrum);
            end
            
            %distant light specific properties
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
        
        %writes the pbrt file corresponding to this object
        function writeFile(obj, fid)
           %call superclass function 
           writeFile@pbrtLightObject(obj, fid);
           
           fprintf(fid,'\t\t"point from" [');
           fprintf(fid,'%f ', obj.from);
           fprintf(fid,']\n');
           
           fprintf(fid,'\t\t"point to" [');
           fprintf(fid,'%f ', obj.to);
           fprintf(fid,']\n');
        end
    end

end