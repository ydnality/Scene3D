classdef ppsfObject < rayObject
    % Create a ray object
    %
    % ray = rayObject(origin,direction,wavelength)
    %
    % Example:
    %   rayObject
    %
    % AL Vistasoft Copyright 2014
    
    properties
%         origin;
%         direction;
%         wavelength;
%         waveIndex;
          pointSourceDepth;          %depth of the originating point source
          pointSourceFieldHeight;    %field height of the originating point source
%           outsideApertureLocation;   %where rays intersect further most lens surface
          apertureSamples;    %where rays intersect the front most aperture
          apertureLocation;     %where rays intersect the actual lens aperture (in the middle usually)
    end
    
    methods
        function obj = ppsfObject(origin, direction, wavelength, pointSourceDepth, pointSourceFieldHeight, apertureSamples)
            %default constructor
            
            % rayObject properties
            if (ieNotDefined('origin')),   obj.origin = [0,0,0];
            else                           obj.origin = origin;
            end
            
            if (ieNotDefined('direction')), obj.direction = [0,0,1]; 
            else                            obj.direction = direction;
            end
            
            if ieNotDefined('wavelength'),  obj.wavelength = 550;
            else                            obj.wavelength = wavelength;
            end
            
            %new properties
            if ieNotDefined('pointSourceDepth'),  obj.pointSourceDepth = 100;
            else                            obj.pointSourceDepth = pointSourceDepth;
            end
            
            if ieNotDefined('pointSourceFieldHeight'),  obj.pointSourceFieldHeight = 0;
            else                            obj.pointSourceFieldHeight = pointSourceFieldHeight;
            end
            
            if ieNotDefined('apertureSamples'),  obj.apertureSamples = 0;
            else                            obj.apertureSamples = apertureSamples;
            end
        end
        
%         function setApertureSamples(obj, apertureSamples)
%             %setApertureLocation(obj, apertureLocation)
%             %
%             %sets the intersection of rays with actual lens aperture
%             %location
%             %TODO: error handling
%             
%            obj.apertureLocation = apertureSamples;
%            return;
%         end

       
        function obj = expandWavelengths(obj, wave, waveIndex)
         %replicates the ray bundle to cover a series of wavelengths   
%             if ieNotDefined('waveIndex'),  waveIndex = 1:length(wave);
%             end
%             subLength = size(obj.origin, 1);
%             obj.origin = repmat(obj.origin, [length(wave) 1]);
%             obj.direction = repmat(obj.direction, [length(wave) 1]);
%             
%             
%             
%             tmp = (wave' * ones(1, subLength))'; %creates a vector representing wavelengths... for example: [400 400 400... 410 410 410... ..... 700]
%             obj.wavelength = tmp(:);
%             
%             %assign the indices of this wavelength expansion - TODO: maybe
%             %make this cleaner somehow...
%             tmp = (waveIndex' * ones(1, subLength))';
%             obj.waveIndex = tmp(:);

            %TODO: see if there's a way to remove this line
            if ieNotDefined('waveIndex'),  waveIndex = 1:length(wave);
            end
            
            obj = expandWavelengths@rayObject(obj, wave, waveIndex);
            
            obj.apertureSamples.X = repmat(obj.apertureSamples.X, [length(wave) 1]);  %added for ppsfObject
            obj.apertureSamples.Y = repmat(obj.apertureSamples.Y, [length(wave) 1]);  %added for ppsfObject
        end
    end
    
end