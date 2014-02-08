% cameraObject contains the camera position, lens, sensor
classdef cameraObject <  handle
    
    properties (SetAccess = private)
        position;
        lens;
        film;    
    end
    methods
        
        %default constructor
        function obj = cameraObject(inPos, inLens, inFilm)
            % Placement of the camera.  Defined by a vector that looks in a certain
            % direction, and then tells you which way is up.
            if (ieNotDefined('inPos'))
                obj.position = ...
                    [4.5 -80 7; % Starting up
                    4.5 -79 7; % Ending up
                    0 0 1];   % Which way is up
            else
                obj.position = inPos;
            end
            
            if (ieNotDefined('inLens'))
                % Example lens
%                 obj.lens = 'idealLens-50mm.pbrt';
                obj.lens = lensPinholeObject();
            else
                validateattributes(inLens, {'lensObject'}, {'nonempty'});
                obj.lens = inLens;
            end
            
            if (ieNotDefined('inFilm'))
                % Typical "sensor"
                obj.film.name = 'image';
                obj.film.xresolution = 200;
                obj.film.yresolution = 200;
            else
                %TODO: verify the film
                obj.film = inFilm;
            end
        end
        
        %TODO: error checking
        function setPosition(obj, inPos)
            obj.position = inPos;
        end
        
        %moves the camera by an offset
        %TODO: error checking
        function moveCamera(obj, offset)
            obj.position = obj.position + cat(1, repmat(offset, [2 1]), [ 0 0 0]);
        end
        
        %sets the lens to inLens
        function setLens(obj, inLens)
            validateattributes(inLens, {'lensObject', 'char'}, {'nonempty'});
            obj.lens = inLens;
        end
 
        %Sets the resolution of the camera to inXres and inYRes
        function setResolution(obj, inXRes, inYRes)
            validateattributes(inXRes, {'numeric'}, {'nonnegative'});
            validateattributes(inYRes, {'numeric'}, {'nonnegative'});
            obj.film.xresolution = inXRes;
            obj.film.yresolution = inYRes;
        end
        
        function returnVal = writeFile(obj, fid)
            %camera position
            fprintf(fid,'\n\nLookAt\n');
            
            if (isfield(obj.position, 'fileName'))
                fprintf(fid,'\n\nInclude "%s"\n', obj.position.fileName);  %TODO: this might need to be fixed later
            else
                fprintf(fid,'\t%f %f %f\n',obj.position');
            end
            
            %lens
            if (ischar(obj.lens))
                fprintf(fid,'\n\nInclude "%s"\n', obj.lens);
            else
                obj.lens.writeFile(fid);
            end
            returnVal = 1;
            
            %film resolution
            fprintf(fid,'\n\nFilm "image"\n');
            fprintf(fid,'\t"integer xresolution" [%i]\n',obj.film.xresolution);
            fprintf(fid,'\t"integer yresolution" [%i]\n',obj.film.yresolution);
            
        end
    end
    
end