% pbrtCameraObject contains the camera position, lens, sensor
classdef pbrtCameraObject < handle
    
    properties (SetAccess = private)
        position;
        transformArray;
        lens; 
        film; 
    end
    
    methods
        
        %default constructor
        function obj = pbrtCameraObject(inPos, inLens, inFilm)
            % Placement of the camera.  Defined by a vector that looks in a certain
            % direction, and then tells you which way is up.
            %
            % Usage of the transformArray and position elements are
            % mutually exlusive.  If TransformArray has elements, then
            % position will be ignored.  If Position is nonempty, then
            % transformArray will be ignored. This is in efforts for
            % backwards compatibility and flexibility.  If both are
            % nonempty, then transformArray has precedence.
            
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
                obj.lens = pbrtLensPinholeObject();
            else
                validateattributes(inLens, {'lensObject'}, {'nonempty'});
                obj.lens = inLens;
            end
            
            if (ieNotDefined('inFilm'))
                % Typical "sensor"
                obj.film.name = 'image';
                obj.film.xresolution = 200;
                obj.film.yresolution = 200;
                obj.film.cropwindow = [0 1 0 1];
            else
                %TODO: verify the film
                obj.film = inFilm;
            end
            
            %initialize the transform array
            obj.transformArray = cell(0,1);
        end
        
        function addTransform(obj, inTrans)
        %addTransform(obj, inTrans)
            if (isa(inTrans, 'pbrtTransformObject'));
                obj.transformArray{end+1} = inTrans;
            else
                error('input transform must be of type pbrtTransformObject');
            end
        end
        
        %TODO: error checking
        function setPosition(obj, inPos)
            %setPosition(obj, inPos)
            obj.position = inPos;
        end
        
        %moves the camera by an offset
        %TODO: error checking
        function moveCamera(obj, offset)
            
            obj.position = obj.position + cat(1, repmat(offset, [2 1]), [ 0 0 0]);
        end
        
        %sets the lens to inLens
        function setLens(obj, inLens)
            validateattributes(inLens, {'pbrtLensObject', 'char'}, {'nonempty'});
            obj.lens = inLens;
        end
 
        %Sets the resolution of the camera to inXres and inYRes
        function setResolution(obj, inXRes, inYRes)
            validateattributes(inXRes, {'numeric'}, {'nonnegative'});
            validateattributes(inYRes, {'numeric'}, {'nonnegative'});
            obj.film.xresolution = inXRes;
            obj.film.yresolution = inYRes;
        end
        
        function setCropWindow(obj, inMinX, inMaxX, inMinY, inMaxY)
        %Sets the crop window of the camera.  All inputs must be between 0
        %and 1.  Crop window is set using ratios of the full film size.
        %
        %inMinX: min bounding window for X
        %inMaxX: max bounding window for X
        %inMinY: min bounding window for Y
        %inMaxy: max bounding window for Y
        
            validateattributes(inMinX, {'numeric'}, {'nonnegative'});
            validateattributes(inMaxX, {'numeric'}, {'nonnegative'});  
            validateattributes(inMinY, {'numeric'}, {'nonnegative'});  
            validateattributes(inMaxY, {'numeric'}, {'nonnegative'});
            
            if (inMinX <= inMaxX && inMinX >= 0 && inMaxX <= 1 && inMinY <= inMaxY && inMinY >=0 && inMaxY <= 1)
                obj.film.cropwindow = [inMinX inMaxX inMinY inMaxY];
            else
               error('Inputs inconsistent! All inputs must be between 0 and 1.  inMinX must be <= inMaxX, and inMinY must be <= inMaxY.'); 
            end
        end
        function returnVal = writeFile(obj, fid)
            
            
            if (isfield(obj.position, 'fileName'))
                fprintf(fid,'\n\nInclude "%s"\n', obj.position.fileName);  %TODO: this might need to be fixed later
            elseif(isempty(obj.transformArray))
                %camera position
                fprintf(fid,'\n\nLookAt\n');
                fprintf(fid,'\t%f %f %f\n',obj.position');
            else
                %we will use transformArray and print those properties
                for i = 1:length(obj.transformArray)
                   obj.transformArray{i}.writeFile(fid); 
                end
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
            fprintf(fid,'\t"float cropwindow" [%f %f %f %f]', obj.film.cropwindow);
            
        end
    end
    
end