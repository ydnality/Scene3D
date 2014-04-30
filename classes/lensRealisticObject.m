classdef lensRealisticObject <  lensObject
    % Create a lens object
    %
    %   lens = lensRealisticObject(elOffset, elRadius, elAperture, elN, aperture, focalLength, center);  % Units are mm
    %
    % Presently we only represent spherical lenses and circular apertures.
    %
    % These are defined by a series of surfaces. We code the offset to each
    % surface radius to the center of the spherical lens.  Positive means
    % sphere center to the right. Aperture parameters (a single number is a
    % diameter in mm). index of refraction (n) for the material to the left
    % of the surface.
    % 
    % pinhole cameras have no aperture and the pinhole lens will inherit
    % this superclas. This will be a superclass that will be inherited by
    % other classes in the future
    %
    % We aim to be consistent with the PBRT lens files, and maybe the Zemax
    % as far possible ?
    %
    % This could become a camera, or we could make a camera object that has
    % a lens and film.
    % 
    % Example:
    %   lensRealisticObject()
    %   lensObject(2, 67, 3, 1, 3, 50, 0)
    %
    % AL Vistasoft Copyright 2014
    
    properties
        elementArray;
        totalOffset;
        numEls;
    end
    
    methods
        
        
        function obj = lensRealisticObject(elOffset, elRadius, elAperture, elN, aperture, focalLength, center, diffractionEnabled, wave)
        %default constructor
            
            % Units are mm
            if (ieNotDefined('elOffset')), elOffset = 0;
            end
            
            % Units are mm
            if (ieNotDefined('elRadius')), elRadius = 10;
            end
            
            % Units are mm
            if (ieNotDefined('elAperture')), elAperture = 10;
            end
            
            % Units are mm
            if (ieNotDefined('elN')), elN = 1;
            end
           
            % Units are mm
            if (ieNotDefined('aperture')), obj.apertureRadius = 3;
            else                           obj.apertureRadius = aperture;
            end            
            
            % Units are mm
            %TODO: error checking
            if (ieNotDefined('center')), obj.centerPosition = [0 0 -1.5];
            else                           obj.centerPosition = center;
            end 
            
            % Units are mm
            %TODO: error checking
            if (ieNotDefined('focalLength')), obj.focalLength = 50;
            else                           obj.focalLength = focalLength;
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
            obj.setElements(elOffset, elRadius, elAperture, elN);
        end
        
        function setElements(obj, elOffset, elRadius, elAperture, elN)
        %sets the lens elements of this realistic lens.  All inputs are
        %vectors and should have equal elements, and be real numbers of
        %type double.

            obj.numEls = length(elOffset); % we must update numEls each time we add a lens element
            
            %error checking
            if (obj.numEls~= length(elRadius) || obj.numEls~= length(elAperture) || obj.numEls ~= size(elN,2))
                error('input vectors must all be of the same lengths');
            end
            
            %if no wavelength dependence of index of refraction specified,
            %(only one column of data, supply one)
            if (size(elN, 1) == 1)
                numWave = length(obj.wave);
                elN = repmat(elN, [numWave 1]);
            end
           
            obj.elementArray = lensElementObject();
            for i = 1:length(elOffset)
                obj.elementArray(i) = lensElementObject(elOffset(i), elRadius(i), elAperture(i), elN(:, i));
            end

            obj.firstApertureRadius = obj.elementArray(1).aperture;
            obj.computeCenters();
            obj.calculateApertureSample(); 
        end
        
        function readLensFile(obj, fullFileName)
        %reads Scene3D format lens files and converts this to our format in
        %the lensRealistic object

            %fid = fopen(fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat'));
            fid = fopen(fullFileName);
            import = textscan(fid, '%s%s%s%s', 'delimiter' , '\t');
            fclose(fid);

            %first find the start of the lens line demarked by #   radius
            firstColumn = import{1};

            %read comment lines
            continu = true;
            i = 1;
            while(continu && i <= length(firstColumn))
                compare = regexp(firstColumn(i), 'radius');
                if(~(isempty(compare{1})))
                    continu = false;
                end
                i = i+1;
            end

            %put data into lens object
            radius = str2double(import{1});
            radius = radius(i:length(firstColumn));
            offset = str2double(import{2});
            offset = offset(i:length(firstColumn));
            %change from pbrt Scene3D format to raytrace Scene3D format
            offset = [0; offset(1:(end-1))];
            N = str2double(import{3});
            N = N(i:length(firstColumn));
            aperture = str2double(import{4})/2;%radius supplied is the radius diameter, so divide it by 2
            aperture = aperture(i:length(firstColumn));
            
            %TODO: eventually calculate this given the lens file
            obj.centerPosition = [0 0 -15.1550];  

            %modify the object and reinitialize
            obj.setElements(offset, radius, aperture, N');
        end
        
        
        function computeTotalOffset(obj)
        %Calculates the total offset of the lens by adding all existing
        %offsets
        %make this private later
            obj.totalOffset  = 0;
            for i = 1:obj.numEls
                obj.totalOffset = obj.totalOffset + obj.elementArray(i).offset;
            end
        end

        
        function computeCenters(obj)
        %computes the spherical centers of each element
        
            obj.computeTotalOffset();
            prevSurfaceZ = -obj.totalOffset;
            for i = 1:length(obj.elementArray)
                obj.elementArray(i).zIntercept = prevSurfaceZ + obj.elementArray(i).offset;  %TODO: these will change later with sets
                obj.elementArray(i).sphereCenter = [0 0 obj.elementArray(i).zIntercept+ obj.elementArray(i).radius];
                prevSurfaceZ = obj.elementArray(i).zIntercept;
            end
        end
        
%         function obj = calculateApertureSample(obj)
%            
%             %loop through aperture positions and uniformly sample the aperture
%             %everything is done in vector form for speed
%             [rectApertureSample.X, rectApertureSample.Y] = meshgrid(linspace(-1, 1, 3),linspace(-1, 1, 3)); %adjust this if needed - this determines the number of samples per light source
%             
%             %assume a circular aperture, and make a mask that is 1 when the pixel
%             %is within a circle of radius 1
%             apertureMask = (rectApertureSample.X.^2 + rectApertureSample.Y.^2) <= 1;
%             scaledApertureSample.X = rectApertureSample.X .* obj.apertureRadius;
%             scaledApertureSample.Y = rectApertureSample.Y .* obj.apertureRadius;
%             
%             %remove cropped sections of aperture
%             croppedApertureSample.X =  scaledApertureSample.X;
%             croppedApertureSample.X(apertureMask == 0) = [];
%             croppedApertureSample.X = croppedApertureSample.X';
%             croppedApertureSample.Y =  scaledApertureSample.Y;
%             croppedApertureSample.Y(apertureMask == 0) = [];
%             croppedApertureSample.Y = croppedApertureSample.Y';
%             
%             obj.apertureSample = croppedApertureSample;
%         end


        function obj =  drawLens(obj)
        %draws the illustration of the lens on a figure - you must declare
        %a new graphwin first!
        
            prevSurfaceZ = -obj.totalOffset;
            prevAperture = 1;

            for lensEl = 1:obj.numEls
                curEl = obj.elementArray(lensEl);

                %illustrations for debug
                if (curEl.radius ~=0)
                    %draw arcs if radius is nonzero
                    zPlot = linspace(curEl.sphereCenter(3) - curEl.radius, curEl.sphereCenter(3) + curEl.radius, 10000);
                    yPlot = sqrt(curEl.radius^2 - (zPlot - curEl.sphereCenter(3)) .^2);
                    yPlotN = -sqrt(curEl.radius^2 - (zPlot - curEl.sphereCenter(3)) .^2);
                    arcZone = 5;
                    %TODO:find a better way to plot the arcs later - this one is prone to potential problem
                    withinRange = and(and((yPlot < curEl.aperture),(zPlot < prevSurfaceZ + curEl.offset + arcZone)), (zPlot > prevSurfaceZ + curEl.offset - arcZone));
                    line(zPlot(withinRange), yPlot(withinRange));
                    line(zPlot(withinRange), yPlotN(withinRange));
                else
                    %draw the main aperture opening if radius = 0
                    
                    %TODO: draw the difference between specified aperture
                    %from file and specified aperture from object instance
                    
                    %right now: take the minimum value
                    
                    curAperture = min(curEl.aperture, obj.apertureRadius);
                    
                    line(curEl.sphereCenter(3) * ones(2,1), [-prevAperture -curAperture]);
                    line(curEl.sphereCenter(3) * ones(2,1), [curAperture prevAperture]);
                end
                
                prevAperture = curEl.aperture;
                prevSurfaceZ = prevSurfaceZ + curEl.offset;                
                
            end
        end


        function obj =  rayTraceThroughLens(obj, rays, debugLines)
            %Ray-trace of the lens
            % The rays are input
            % The rays are updated after refraction

            % The order is from furthest from film to film, which is also
            % how the rays pass through the optics.
            
            %
            
            passedCenterAperture = false;  %true if rays are traced through lens aperture
            if (ieNotDefined('debugLines'))
                debugLines = false;
            end
            
            prevSurfaceZ = -obj.totalOffset;
            prevN = ones(length(rays.origin), 1);
            
            obj.drawLens();
                         
            for lensEl = 1:obj.numEls
                curEl = obj.elementArray(lensEl);
                curAperture = curEl.aperture;
                
                %  ----vectorized               
                
                %ray trace through a single element
                %calculate intersection with lens element -
                if (curEl.radius ~= 0) %only do this for actual spherical elements, 
                    repCenter = repmat(curEl.sphereCenter, [length(rays.origin) 1]);
                    repRadius = repmat(curEl.radius, [length(rays.origin) 1]);
                    radicand = dot(rays.direction, rays.origin - repCenter, 2).^2 - ...
                        ( dot(rays.origin - repCenter, rays.origin -repCenter, 2)) + repRadius.^2;
                    if (curEl.radius < 0)
                        intersectT = (-dot(rays.direction, rays.origin - repCenter, 2) + sqrt(radicand));
                    else
                        intersectT = (-dot(rays.direction, rays.origin - repCenter, 2) - sqrt(radicand));
                    end

                    %make sure that T is > 0
                    %                         if (intersectT < 0)
                    %                             disp('Warning: intersectT less than 0.  Something went wrong here...');
                    %                         end

                    repIntersectT = repmat(intersectT, [1 3]);
                    intersectPosition = rays.origin + repIntersectT .* rays.direction;

                    %illustrations for debugging
                    if (debugLines)
                        xCoordVector = [rays.origin(:,3) intersectPosition(:,3) NaN([length(rays.origin) 1])]';
                        yCoordVector = [rays.origin(:,2) intersectPosition(:,2) NaN([length(rays.origin) 1])]';
                        xCoordVector = real(xCoordVector(:));
                        yCoordVector = real(yCoordVector(:));
                        line(xCoordVector,  yCoordVector ,'Color','b','LineWidth',1);
                    end
                        
                    %                         if (isnan(intersectPosition))
                    %                             disp('nan value');
                    %                         end
                else       
                    %plane intersection with lens aperture - TODO: maybe put
                    %in function?
                    intersectZ = repmat(curEl.sphereCenter(3), [length(rays.origin) 1]); %assumes that aperture is perfectly perpendicular with optical axis
                    intersectT = (intersectZ - rays.origin(:, 3))./rays.direction(:, 3);
                    repIntersectT = repmat(intersectT, [1 3]);
                    intersectPosition = rays.origin + rays.direction .* repIntersectT;
                    curAperture = min(curEl.aperture, obj.apertureRadius);
                    
                    %added for ppsfObject apertureTracking
                    if(isa(rays, 'ppsfObject'))
                        rays.setApertureLocation(intersectPosition);
                        passedCenterAperture = true;
                    end
                    %                         if (isnan(intersectPosition))
                    %                             disp('nan value');
                    %                         end
                end
 
                % remove rays that land outside of the aperture
                %TODO: consider making these set functions later
                outsideAperture = intersectPosition(:, 1).^2 + intersectPosition(:, 2).^2 > curAperture^2;
                rays.origin(outsideAperture, : ) = [];
                rays.direction(outsideAperture, : ) = [];
                rays.wavelength(outsideAperture) = [];
                rays.waveIndex(outsideAperture) = [];
                intersectPosition(outsideAperture, :) = [];
                prevN(outsideAperture) = [];
                
                %special case with ppsfObjects
                if(isa(rays,'ppsfObject') && passedCenterAperture)
                    rays.apertureLocation(outsideAperture, :) = [];   
                end
                
                % snell's law
                if(curEl.radius ~= 0)
                    %in bounds case - perform vector snell's law
                    repCenter = repmat(curEl.sphereCenter, [length(rays.origin) 1]);
                    normalVec = intersectPosition - repCenter;  %does the polarity of this vector matter? YES
                    normalVec = normalVec./repmat(sqrt(sum(normalVec.*normalVec, 2)),[1 3]); %normalizes each row 

                    if (curEl.radius < 0)  %which is the correct sign convention? This is correct
                        normalVec = -normalVec;
                    end

                    %modify the index of refraction depending on wavelength
                    %TODO: have this be one of the input parameters (N vs. wavelength)
%                     if (curEl.n ~= 1)
%                         curN = (rays.wavelength - 550) * -.04/(300) + curEl.n;
%                     else
%                         curN = ones(length(rays.wavelength), 1);
%                     end
                    
%                     waveIndex = find([obj.wave' obj.nWave']== rays.wavelength);  %this doesn't work yet
                    %this was supposed to conver wavelength to waveIndex -
                    %but I coudln't find a way to vectorize it, so instead
                    %we precompute it
                    
                    curN = curEl.n(rays.waveIndex);
                    
                    
%                     curN = ones(length(rays.wavelength), 1) * curEl.n;
                    ratio = prevN./curN;    %snell's law index of refraction

                    %Vector form of Snell's Law
                    c = -dot(normalVec, rays.direction, 2);
                    repRatio = repmat(ratio, [1 3]);
                    newVec = repRatio .* rays.direction + repmat((ratio.*c -sqrt(1 - ratio.^2 .* (1 - c.^2))), [1 3])  .* normalVec;
                    newVec = newVec./repmat(sqrt(sum(newVec.*newVec, 2)), [1 3]); %normalizes each row

                    %update the direction of the ray
                    rays.origin = intersectPosition;
                    rays.direction = newVec;
                    prevN = curN;  %note: curN won't change if the aperture is the overall lens aperture
                end

                %HURB diffraction
                if (obj.diffractionEnabled)
                    obj.rayTraceHURB(rays, intersectPosition, curEl.aperture);
                end


                %iterate previous z 
                prevSurfaceZ = prevSurfaceZ + curEl.offset;
            end
        end     
        
        
        
        
        function obj =  rayTraceThroughLensNonvectorized(obj, rays)
        %performs ray-trace of the lens, given an input bundle or rays
        %outputs the rays that have been refracted by the lens
        %This function is here mostly for reference and debugging purposes
        %- it is now replaced by rayTraceThroughLens which is the
        %vectorized version
        %the order of the lens elements has since been reversed - here is
        %an example input
        % % offset = [1.5 1.5 0];
        % radius = [-67 0 67];
        % aperture = [5 1 4];
        %n = [ 1 0 1.67];
        
            %initialize newRays to be the old ray.  We will update it later.
%             newRays = rays;
            
            prevSurfaceZ = -obj.totalOffset;
            prevN = 1;
%             
%             prevN = ones(length(rays.origin), 1);  %assume that we start off in air
%             
            for lensEl = obj.numEls:-1:1
                curEl = obj.elementArray(lensEl);
%                 curEl.center = [0 0 prevSurfaceZ + curEl.offset + curEl.radius];
                
                %illustrations for debug
                zPlot = linspace(curEl.sphereCenter(3) - curEl.radius, curEl.sphereCenter(3) + curEl.radius, 10000);
                yPlot = sqrt(curEl.radius^2 - (zPlot - curEl.sphereCenter(3)) .^2);
                yPlotN = -sqrt(curEl.radius^2 - (zPlot - curEl.sphereCenter(3)) .^2);
                arcZone = 5;
                %TODO:find a better way to plot the arcs later - this one is prone to potential problem
                withinRange = and(and((yPlot < curEl.aperture),(zPlot < prevSurfaceZ + curEl.offset + arcZone)), (zPlot > prevSurfaceZ + curEl.offset - arcZone));
                line(zPlot(withinRange), yPlot(withinRange));
                line(zPlot(withinRange), yPlotN(withinRange));

%                   -nonvectorized loop ----   
                
                
%                 newRays = rayObject(rays.origin, rays.direction, rays.wavelength);
                i = 1;
                while(i <= length(rays.origin))
                    %get the current ray
                    
                    ray = rayObject(rays.origin(i,:), rays.direction(i,:), rays.wavelength(i));
                    
                    if (isnan(ray.direction))
                        disp('nan value');
                     end
%                     ray.direction = rays.direction(i,:);   %TODO: replace with real ray object
%                     ray.origin = rays.origin(i,:);
%                     ray.wavelength = rays.wavelength(i);
                    
                    %calculate intersection with lens element -
                    %TODO: maybe put in function?
                    if (curEl.radius ~= 0) %only do this for actual spherical elements, 
                        radicand = dot(ray.direction, ray.origin - curEl.sphereCenter)^2 - ...
                            ( dot(ray.origin -curEl.sphereCenter, ray.origin -curEl.sphereCenter)) + curEl.radius^2;
                        if (curEl.radius < 0)
                            intersectT = (-dot(ray.direction, ray.origin - curEl.sphereCenter) + sqrt(radicand));
                        else
                            intersectT = (-dot(ray.direction, ray.origin - curEl.sphereCenter) - sqrt(radicand));
                        end
                        
                        %make sure that T is > 0
                        if (intersectT < 0)
                            disp('Warning: intersectT less than 0.  Something went wrong here...');
                        end
                        
                        intersectPosition = ray.origin + intersectT * ray.direction;
                        
                        %illustrations for debugging
                        %                     lensIllustration(max(round(intersectPosition(2) * 100 + 150),1), max(-round(intersectPosition(3) * 1000), 1)) = 1;  %show a lens illustration
                        hold on;
                        %TODO: potential problem with non-real answers here
                        line(real([ray.origin(3) intersectPosition(3) ]), real([ray.origin(2) intersectPosition(2)]) ,'Color','b','LineWidth',1);
                        
                        if (isnan(intersectPosition))
                            disp('nan value');
                        end
                    else       
                        %plane intersection with lens aperture - TODO: maybe put
                        %in function?
                        intersectZ = curEl.sphereCenter(3); %assumes that aperture is perfectly perpendicular with optical axis
                        intersectT = (intersectZ - ray.origin(:, 3))./ray.direction(:, 3);
                        intersectPosition = ray.origin + ray.direction * intersectT;
                        
                        if (isnan(intersectPosition))
                            disp('nan value');
                        end
                    end
                    
                    
%                     check for bounds within aperture, 
%                     if(intersectPosition(1)^2 + intersectPosition(2)^2 > curEl.aperture^2)
% %                         ray.origin = [];  %not in bounds - remove
% %                         ray.direction = [];
% %                         ray.wavelength = [];
%                         %not in bounds - remove
%                         rays.origin(i, : ) = [];
%                         rays.direction(i, : ) = [];
%                         rays.wavelength(i) = [];
%                     else
                        if(curEl.radius ~= 0) %only perform Snell's law if not the lens aperture
                        
                            %in bounds case - perform vector snell's law
                            normalVec = intersectPosition - curEl.sphereCenter;  %does the polarity of this vector matter? YES
%                             normalVec = normalVec./norm(normalVec);
                            normalVec = normalVec./sqrt(sum(normalVec .* normalVec));
                            
                            if (curEl.radius < 0)  %which is the correct sign convention? This is correct
                                normalVec = -normalVec;
                            end
                            
                            %modify the index of refraction depending on wavelength
                            %TODO: have this be one of the input parameters (N vs. wavelength)
                            if (curEl.n ~= 1)
                                curN = (ray.wavelength - 550) * -.04/(300) + curEl.n;
                            else
                                curN = 1;
                            end
                            
                            %this step is in here to correct for the
                            %previous element N depending on wavelength
                            if (prevN ~= 1)
                                prevNAdjusted = (ray.wavelength - 550) * -.04/(300) + prevN;
                            else
                                prevNAdjusted = 1;
                            end
                            ratio = prevNAdjusted/curN;    %snell's law index of refraction
                            
                            %Vector form of Snell's Law
                            c = -dot(normalVec, ray.direction);
                            newVec = ratio *ray.direction + (ratio*c -sqrt(1 - ratio^2 * (1 - c^2)))  * normalVec;
%                             newVec = newVec./norm(newVec); %normalize
                             newVec = newVec./sqrt(sum(newVec .* newVec)); %normalize
                            
                            %update the direction of the ray
                            ray.origin = intersectPosition;
                            ray.direction = newVec;
                            
                            if (isnan(ray.direction))
                              disp('nan value');  
                            end
                        end
                        
                        
                        
                        
                        %                     tempRay = rayObject(rays.origin(i, : ) , rays.direction(i, : ), rays.wavelength(i));
                        % diffraction HURB calculation
                        if (obj.diffractionEnabled)
                            obj.rayTraceHURB(ray, intersectPosition, curEl.aperture);
                            if (isnan(ray.direction))
                              disp('nan value');  
                            end
                        end
                        
                        rays.origin(i, : ) = ray.origin;
                        rays.direction(i, : ) = ray.direction;
                        rays.wavelength(i) = ray.wavelength;
                        i = i +1;
                    end
       
%                 end
                 % -nonvectorized loop ----   
                
                %remove dead rays
%                 rays.origin(rays.origin == -999) = [];
%                 rays.direction(rays.direction == -999) = [];
%                 rays.wavelength(rays.wavelength == -999) = [];
                
                %iterate index of refraction and previous z 
                if(curEl.radius ~= 0)
                    prevN = curEl.n;
                end
                prevSurfaceZ = prevSurfaceZ + curEl.offset;
            end
        end     
        
    end
    
end