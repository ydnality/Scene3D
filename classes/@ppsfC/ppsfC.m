classdef ppsfC < rayC
    % Create a plenoptic point spread function object
    %
    %   ppsf = ppsfC(origin,direction,wavelength)
    %
    % This is related to ray objects, must explain more here.
    %
    % Example:
    %   ppsf = ppsfC;
    %
    % TODO: in order to save memory, we anticipate that all the
    % intersection data members will contain a Z intercept separately, and
    % also an X,Y intercept.  This is because all the intercepts will
    % contain the identical Z intercept.
    %
    % AL Vistasoft Copyright 2014
    
    properties
          %location of the originating point source scene is -Z direction
          pointSourceLocation = [ 0 0 -100];
          
          %apertureSamples - where rays intersect the front most aperture
          aEntranceInt = struct('XY', [0 0], 'Z', 0);  
          
          %apertureLocation - where rays intersect the actual lens aperture (in the middle usually)
          aMiddleInt = struct('XY', [0 0], 'Z', 0);   
          
          % Exit positions
          aExitInt = struct('XY', [0 0], 'Z', 0);
          
          %exit direction of light field
          aExitDir = 0;   
    end
    
    methods
        function obj = ppsfC(varargin)
            % Constructor
            
            % Set up defaults.
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
        end
               
        function obj = expandWavelengths(obj, wave, waveIndex)
         %replicates the ray bundle to cover a series of wavelengths   

            %TODO: see if there's a way to remove this line
            if ieNotDefined('waveIndex'),  waveIndex = 1:length(wave);
            end
            
            obj = expandWavelengths@rayC(obj, wave, waveIndex);
                       
            obj.aEntranceInt.XY = repmat(obj.aEntranceInt.XY, [length(wave) 1]);  %added for ppsfC
        end
        
        
        function obj = projectOnPlane(obj, planeLocation)
            %projectOnPlane(obj, planeLocation)
            %planeLocation: the z coordinate of a plane that is parallel to the
            %x-y plane.
            %intersects the rays with a plane at a specified location.
            %This is meant to make it easier to analyze the lightfield function.
            
            %if we are projecting onto the z = 0 plane, assign the exit
            %pupil intersection data member
            projectOnPlane@rayC(obj, planeLocation);
            
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