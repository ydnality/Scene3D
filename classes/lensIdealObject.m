classdef lensIdealObject <  lensObject
    % Create a lens object
    %
    %   lens = lensIdealObject(aperture,focalLength);  % Units are mm
    %
    % An ideal thin lens is defined as one that obey's the thin lens equation
    % (1/s1 + 1/s2 = 1/f).  Therefore, Snell's Law is not used in this
    % object at all.  
    % Instead, the direction that the rays bend are determined by the
    % thin lens equation.  At each point source, a "center ray" is shot at the
    % center of the lens.  The intersection of this ray and the focal-plane as
    % defined by the thin lens equation determines the point of focus of the
    % lens.  All other rays that are shot at the edge of the aperture will then
    % intersect this ideal focal point.
    % All units are in mm.
    % 
    % Example:
    %   lensIdealObject()
    %   lensIdealObject(3, 50)
    %
    % AL Vistasoft Copyright 2014
    
    properties
        
    end
    
    methods
        
        %default constructor - TODO:  figure out some inheritance issues
        %for efficiency
        function obj = lensIdealObject(aperture, focalLength, center, diffractionEnabled, wave)
            
            % Units are mm
            if (ieNotDefined('aperture')), obj.apertureRadius = 3;
            else                           obj.apertureRadius = aperture;
            end            
            
            % Units are mm
            if (ieNotDefined('focalLength')), obj.focalLength = 50;
            else                           obj.focalLength = focalLength;
            end          
            
             % Units are mm
            %TODO: error checking
            if (ieNotDefined('center')), obj.centerPosition = [0 0 0];
            else                           obj.centerPosition = center;
            end 
            
            %TODO: error checking
            if (ieNotDefined('diffractionEnabled')), obj.diffractionEnabled = false;
            else                           obj.diffractionEnabled = diffractionEnabled;
            end 
            
            %TODO: error checking
            if (ieNotDefined('wave')), obj.wave = 400:10:700;
            else                           obj.wave = wave;
            end          
            
            obj.nWave = 1:length(obj.wave);
            obj.firstApertureRadius = obj.apertureRadius;
            obj.calculateApertureSample();
        end

        
        
        function obj =  rayTraceThroughLens(obj, rays, curPointSource)
        %traces rays through the lens
            
            % --- center ray calculation-------------
            %trace ray from point source to lens center, to image.  This helps
            %determine the point of focus
            centerRay.origin = curPointSource;
            centerRay.direction = obj.centerPosition - centerRay.origin;
            centerRay.direction = centerRay.direction./norm(centerRay.direction);
            
            %calculate the in-focus plane using thin lens equation
            inFocusDistance = 1/(1/obj.focalLength - -1/curPointSource(3));
            
            %calculates the in-focus position.  The in-focus position is the
            %intersection of the in-focus plane and the center-ray
            inFocusT = (inFocusDistance - centerRay.origin(3))/centerRay.direction(3);
            inFocusPosition = centerRay.origin + inFocusT .* centerRay.direction;
            % --------center ray calculation -------
            
            
            %---------------- lens refraction code  -----
            %when intersecting ideal lens, change the direction to intersect the
            %inFocusPosition, and update the origin
            lensIntersectT = (obj.centerPosition(3) - rays.origin(:,3))./ rays.direction(:,3);
            lensIntersectPosition = rays.origin +  repmat(lensIntersectT, [1 3]) .* rays.direction;
            
            %debug visualization
            for i = 1:size(rays.direction,1)
                hold on;
                line([rays.origin(i,3) lensIntersectPosition(i,3) ], [rays.origin(i,2) lensIntersectPosition(i,2)] ,'Color','b','LineWidth',1);
            end
            
            %calculate new direction
            %             newRays = rayObject(); % added
            rays.origin = lensIntersectPosition;
            rays.direction = repmat(inFocusPosition , [size(rays.origin,1) 1 ]) - rays.origin;
            rays.wavelength = rays.wavelength;
            
            % diffraction HURB calculation
            if (obj.diffractionEnabled)
                obj.rayTraceHURB(rays, lensIntersectPosition, obj.apertureRadius);
            end
            %---------------- lens refraction code  -----
            
        end     
    end
    
end