classdef ppsfObject < rayObject
    % Create a ray object
    %
    % ppsf = ppsfObject(origin,direction,wavelength)
    %
    % Example:
    %   ppsf
    %
    % TODO: in order to save memory, we anticipate that all the
    % intersection data members will contain a Z intercept separately, and
    % also an X,Y intercept.  This is because all the intercepts will
    % contain the identical Z intercept.
    %
    % AL Vistasoft Copyright 2014
    
    properties
%         origin;
%         direction;
%         wavelength;
%         waveIndex;
          pointSourceLocation = [ 0 0 -100];          %location of the originating point source scene is -Z direction
%           pointSourceFieldHeight = 0;    %field height of the originating point source
%           outsideApertureLocation;   %where rays intersect further most lens surface

          aEntranceInt = struct('XY', [0 0], 'Z', 0);  %apertureSamples;    %where rays intersect the front most aperture
          aMiddleInt = struct('XY', [0 0], 'Z', 0);   %apertureLocation;     %where rays intersect the actual lens aperture (in the middle usually)
          aExitInt = struct('XY', [0 0], 'Z', 0);
          aExitDir = 0;   %exit direction of light field
    end
    
    methods
        function obj = ppsfObject(varargin)
            
%             (origin, direction, wavelength, pointSourceDepth, pointSourceFieldHeight, aEntranceInt)
            %default constructor
            
           for ii=1:2:length(varargin)
                p = ieParamFormat(varargin{ii});
                switch p
                    case 'origin'
                        obj.origin = varargin{ii+1};
                    case 'direction'
                        obj.direction = varargin{ii+1};  %must be a 2 element vector
                    case 'wavelength'
                        obj.wavelength = varargin{ii+1};
                    case 'pointsourcedepth'
                         obj.pointSourceDepth = varargin{ii+1};
                    case 'pointsourcefieldheight'
                         obj.pointSourceFieldHeight = varargin{ii+1};
                    case 'aentranceint'
                        obj.aEntranceInt = varargin{ii+1};
                    otherwise
                        error('Unknown parameter %s\n',varargin{ii});
                end
           end 
            
                        
            
            
            % rayObject properties
%             if (ieNotDefined('origin')),   obj.origin = [0,0,0];
%             else                           obj.origin = origin;
%             end
%             
%             if (ieNotDefined('direction')), obj.direction = [0,0,1]; 
%             else                            obj.direction = direction;
%             end
%             
%             if ieNotDefined('wavelength'),  obj.wavelength = 550;
%             else                            obj.wavelength = wavelength;
%             end
%             
%             %new properties
%             if ieNotDefined('pointSourceDepth'),  obj.pointSourceDepth = 100;
%             else                            obj.pointSourceDepth = pointSourceDepth;
%             end
%             
%             if ieNotDefined('pointSourceFieldHeight'),  obj.pointSourceFieldHeight = 0;
%             else                            obj.pointSourceFieldHeight = pointSourceFieldHeight;
%             end
%             
%             if ieNotDefined('aEntranceInt'),  obj.aEntranceInt = 0;
%             else                            obj.aEntranceInt = aEntranceInt;
%             end
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
            
%             obj.aEntranceInt.X = repmat(obj.aEntranceInt.X, [length(wave) 1]);  %added for ppsfObject
%             obj.aEntranceInt.Y = repmat(obj.aEntranceInt.Y, [length(wave) 1]);  %added for ppsfObject
            
            obj.aEntranceInt.XY = repmat(obj.aEntranceInt.XY, [length(wave) 1]);  %added for ppsfObject
%             obj.aEntranceInt.Y = repmat(obj.aEntranceInt.Y, [length(wave) 1]);  %added for ppsfObject
        end
        
        
        function obj = projectOnPlane(obj, planeLocation)
            %projectOnPlane(obj, planeLocation)
            %planeLocation: the z coordinate of a plane that is parallel to the
            %x-y plane.
            %intersects the rays with a plane at a specified location.
            %This is meant to make it easier to analyze the lightfield function.
            
            %if we are projecting onto the z = 0 plane, assign the exit
            %pupil intersection data member
            projectOnPlane@rayObject(obj, planeLocation);
            
            %special case for z = 0 plane - then we modify the ppsf
            %structure for exit aperture intersection position
            if (planeLocation == 0)
                obj.aExitInt.XY = 0;
                obj.aExitInt.XY = obj.origin(:,1:2);    %take the XY coords
                obj.aExitInt.Z = 0;  %only 1 scalar here
                obj.aExitDir = obj.direction;
            end
            
        end
    end
    
end