classdef lensC <  handle
    % Create a multiple element lens object
    %
    %   lens = lensC(elOffset, elRadius, elAperture, elN, aperture, focalLength, center);  
    %
    % Distance units, when not specified, are millimeters.
    %
    % Presently we represent multi-element lenses as a set of spherical
    % lenses and circular apertures. The multiple elements are defined by a
    % series of surfaces with a curvature, position, and index of
    % refraction (as a function of wavelength).
    %
    % We code the offset to each surface radius to the center of the
    % spherical lens.  Positive means sphere center to the right. Aperture
    % parameters (a single number is a diameter in mm). index of refraction
    % (n) for the material to the left of the surface.
    %
    % Lens files and surface arrays are specified from the front most
    % element (closest to the scene), to the back-most element (closest to
    % the sensor)
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
    % Todo: consider a set function for better error handling.
    %
    % Example:
    %   lensC()
    %
    % AL Vistasoft Copyright 2014
    
    properties
        name = 'default';
        type = 'multi element lens';            % We might be removing this soon (why? - BW)
        surfaceArray = surfaceC();     % Set of spherical surfaces and apertures
        diffractionEnabled = false;% Not implemented yet
        wave = 400:50:700;         % nm
        focalLength = 50;          % mm, focal length of multi-element lens
        apertureMiddleD = 1;       % mm, diameter of the middle aperture
        apertureSample = [11 11];  % Number of spatial samples in the aperture.  Use odd number
        centerZ = 0;    %theoretical center of lens (length-wise) in the z coordinate
        
        centerRay = [];   %for use for ideal Lens
        inFocusPosition = [0 0 0];  
    end
    
    methods (Access = public)
        %Multiple element lens constructor
        %TODO: error handling
        function obj = lensC(varargin)
            %  surfaceList = lensReadFile(fName);
            %  lensME = lensC('surfaceArray',surfList);
            
            for ii=1:2:length(varargin)
                p = ieParamFormat(varargin{ii});
                switch p
                    case 'name'
                        obj.name = varargin{ii+1};
                    case 'type'
                        obj.type = varargin{ii+1};
                    case 'surfacearray'
                        obj.surfaceArray = varargin{ii+1};
                    case 'aperturesample'
                        obj.apertureSample = varargin{ii+1};  %must be a 2 element vector
                    case 'aperturemiddled'
                        obj.apertureMiddleD = varargin{ii+1};
                    case 'focallength'
                        obj.focalLength = varargin{ii+1};
                    case 'diffractionenabled'
                        obj.diffractionEnabled = varargin{ii+1};
                    case 'wave'
                        obj.wave = varargin{ii+1};
                    case 'filename'
                        obj.fileRead(varargin{ii+1});
                    otherwise
                        error('Unknown parameter %s\n',varargin{ii});
                end
            end
        end
        
        % Get properties
        function res = get(obj,pName,varargin)
            % Get various derived lens properties though this call
            pName = ieParamFormat(pName);
            switch pName
                case 'name'
                    res = obj.name;
                case 'type'
                    res = obj.type;
                case {'nsurfaces','numels'}
                    % Should be nsurfaces
                    res = length(obj.surfaceArray);
                case 'totaloffset'
                    res = -(obj.surfaceArray(1).sCenter(3) - obj.surfaceArray(1).sRadius);
                case 'surfacearray'
                    res = obj.surfaceArray;
                case 'sradius'
                    % spherical radius of curvature of this surface.
                    % lens.get('sradius',whichSurface)
                    if isempty(varargin), this = 1;
                    else this = varargin{1};
                    end
                    res = obj.surfaceArray(this).sRadius; 
               case 'sdiameter'
                    % Aperture diameter of this surface.
                    % lens.get('sradius',whichSurface)
                    if isempty(varargin), this = 1;
                    else this = varargin{1};
                    end
                    res = obj.surfaceArray(this).apertureD;
                otherwise
                    error('Unknown parameter %s\n',pName);
            end
        end
        
        %% Set
        % This could go away some day.  But for now, the wavelength set is
        % ridiculous because there are so many copies of wave.  So, we
        % should set it here rather than addressing every surface element.
        function set(obj,pName,val,varargin)
            pName = ieParamFormat(pName);
            switch pName
                case 'wave'
                    % The wavelength is annoying.
                    obj.wave = val;
                    nSurfaces = obj.get('n surfaces');
                    for ii=1:nSurfaces
                        obj.surfaceArray(ii).set('wave', val);
                    end
                case 'nall'
                    % Set the index of refraction to all the surfaces
                    nSurfaces = obj.get('n surfaces');
                    for ii=1:nSurfaces
                        obj.surfaceArray(ii).n = val;
                    end    
                otherwise
                    error('Unknown parameter %s\n',pName);
            end
            
        end
        
    end

    
    % Not sure why these are private. (BW).
    methods (Access = private)
        
        function apertureMask = apertureMask(obj)
            % Identify the grid points on the circular aperture.
            %
            
            % Here is the full sampling grid for the resolution and
            % aperture radius
            aGrid = obj.fullGrid;
            
            % We assume a circular aperture. This mask is 1 for the sample
            % points within a circle of the aperture radius
            firstApertureRadius = obj.surfaceArray(1).apertureD/2;
            apertureMask = (aGrid.X.^2 + aGrid.Y.^2) <= firstApertureRadius^2;
            % vcNewGraphWin;  mesh(double(apertureMask))
            
        end
        
        function aGrid = fullGrid(obj,randJitter)
            % Build the full sampling grid, possibly adding a little jitter
            % to avoid aliasing artifacts
            
            if (ieNotDefined('randJitter')), randJitter = false; end
            
            % First make the rectangular samples.
            firstApertureRadius = obj.surfaceArray(1).apertureD/2;
            xSamples = linspace(-firstApertureRadius, firstApertureRadius, obj.apertureSample(1));
            ySamples = linspace(-firstApertureRadius, firstApertureRadius, obj.apertureSample(2));
            [X, Y] = meshgrid(xSamples,ySamples);
            
            %Add a random jitter.  Uniform distribution.  Scales to plus or
            %minus half the sample width
            if(randJitter)
                X = X + (rand(size(X)) - .5) * (xSamples(2) - xSamples(1));
                Y = Y + (rand(size(Y)) - .5) * (ySamples(2) - ySamples(1));
            end
            aGrid.X = X; aGrid.Y = Y;
            
        end
        
        function aGrid = apertureGrid(obj,randJitter)
            % Find the full grid, mask it and return only the (X,Y)
            % positions inside the masked region.  This is the usual set of
            % positions that we use for calculating light fields.
            if ieNotDefined('randJitter'), randJitter = false; end
            
            aGrid = fullGrid(obj,randJitter);
            aMask = apertureMask(obj);
            aGrid.X = aGrid.X(aMask);
            aGrid.Y = aGrid.Y(aMask);
            
            %debug check
            % figure; plot(aGrid.X, aGrid.Y, 'o');
        end
        
        function centers = centersCompute(obj, sOffset, sRadius)
            % Computes the center position of each spherical lens element
            % from the array of element offsets and radii. 
            %
            % N.B.  obj is not used here.  Maybe this shouldn't be in
            % the class definition?
            nSurfaces = length(sOffset);
            
            % The lens file includes the offsets between each surface and its
            % previous one (millimeters).  The location of the first surface is
            % at negative of the total offset.
            zIntercept = zeros(nSurfaces,1);
            zIntercept(1) = -sum(sOffset);
            
            % Cumulte the offsets to find the intercept of the following
            % surface
            for ii=2:length(sOffset)
                zIntercept(ii) = zIntercept(ii-1) + sOffset(ii);
            end
            
            % Centers of the spherical surfaces are stored in this matrix
            % (nSurfaces x 3)
            z = zeros(nSurfaces,1);
            centers = [z z zIntercept + sRadius];  
        end
        
        function elementsSet(obj, sOffset, sRadius, sAperture, sN)
            %Copies the lens elements into the surfaceArray slot of the
            %multielement lens object.
            %
            %The input vectors and should have equal length, and be of type
            %double.
            %
            % Should this just be called from lens.set('surface array')?
            
            
            % Check argument size
            if (length(sOffset)     ~= length(sRadius) || ...
                    length(sOffset) ~= length(sAperture) || ...
                    length(sOffset) ~= size(sN,1))
                error('Input vectors must be of equal length');
            end
            
            % If no wavelength dependence of index of refraction specified,
            % we assume constant across all measurement wavelengths
            % sN is nSurfaces x nWave
            if (size(sN,2) == 1)
                sN = repmat(sN, [1 length(obj.wave)]);
            end

            
            %create array of surfaces
            obj.surfaceArray = surfaceC();
            
            %compute surface array centers
            centers = obj.centersCompute(sOffset, sRadius);
            
            for i = 1:length(sOffset)
                obj.surfaceArray(i) = ...
                    surfaceC('sCenter', centers(i, :), ...
                    'sRadius', sRadius(i), ...
                    'apertureD', sAperture(i), 'n', sN(i, :));
            end
            
        end
        
        function obj = rtRealisticThroughLens(obj, rays, nLines)
            % lens.rtRealisticThroughLens(rays,nLines)
            % 
            %   ray structure 
            %   nLines specifies rendering options
            %     This can be a structure with the fields
            %       .spacing ('uniform' or 'random')
            %       .numLines (how many lines)
            %     This can be a number
            %       <= 0 means don't draw the lines; any other positive
            %     number describes how many to draw.
            %
            % The initial rays are generated by a call to
            % lens.rtSourceToEntrance, which takes a point input and
            % generates the ray positions and angles at the entrance plane.
            %
            % On return, these rays are changed to be the position and
            % direction of the rays at the exit aperture.
            %
            % This should be handled, in the end, by a linear transform
            % following our analysis of the lens ABCD matrices and
            % lightfields.
            %
            % See also psfCameraC.estimatePSF
            %   In that routine, one can ask that the lines be extended to
            %   the film plane by setting a flag.  The flag adds a final
            %   surface in the film plane to the surfaceArray 
            %
            % TODO:  Simplify this code and add comments.
            %        Especially, use the Wigner/Light field ideas
            
            
            % Ray trace calculation starts here
            %
            % The order is from furthest from film to film, which is also
            % how the rays pass through the optics.
            
            lWidth = 0.1; lColor = [0 0.5 1]; lStyle = '-';
            % passedCenterAperture = false;  %true if rays are traced through lens aperture
             
            nRays = rays.get('n rays');
            
            % prevSurfaceZ = -obj.get('totalOffset');
            prevN = ones(nRays, 1);
            
            % For each surface element (lenses and apertures).
            nSurfaces = obj.get('numels');
            
            for lensEl = 1:nSurfaces
                
                % Get the surface data
                curEl = obj.surfaceArray(lensEl);
                curAperture = curEl.apertureD/2;
 
                % Calculate ray intersection position with lens element or
                % aperture. In the case of a 0 curvature, the direction
                % does not change.
                %
                % This uses the vector form of Snell's Law:
                % http://en.wikipedia.org/wiki/Snell's_law
                
                if (curEl.sRadius ~= 0)
                    
                    % Spherical element
                    repCenter = repmat(curEl.sCenter, [nRays 1]);
                    repRadius = repmat(curEl.sRadius, [nRays 1]);
                    
                    % Radicand from vector form of Snell's Law
                    radicand = dot(rays.direction, rays.origin - repCenter, 2).^2 - ...
                        ( dot(rays.origin - repCenter, rays.origin -repCenter, 2)) + repRadius.^2;
                    
                    % Calculate something about the ray angle with respect
                    % to the current surface.  AL to figure this one out
                    % and put in a book reference.
                    if (curEl.sRadius < 0)
                        intersectT = (-dot(rays.direction, rays.origin - repCenter, 2) + sqrt(radicand));
                    else
                        intersectT = (-dot(rays.direction, rays.origin - repCenter, 2) - sqrt(radicand));
                    end
                    
                    %make sure that intersectT is > 0
                    if (min(intersectT(:)) < 0)
                        fprintf('intersectT less than 0 for lens %i',lensEl);
                    end
                    
                    repIntersectT = repmat(intersectT, [1 3]);
                    intersectPosition = rays.origin + repIntersectT .* rays.direction;
                    
                    % Let's try to make the lines 3D.  That would be cool!
                    
                    %nLines can be either a number or a structure
                    if (isfield(nLines, 'spacing') && isfield(nLines,'numLines'))
                        % Structure case.
                        %   .spacing is either 'uniform' or 'random'.  
                        %   .numLines is a positive integer of rays to draw
                        if (nLines.numLines > 0)
                            if lensEl ==1
                                rays.plotHandle = vcNewGraphWin;
                                obj.draw();   % Draw the lens
                                if (strcmp(nLines.spacing, 'uniform'))
                                    samps = round(linspace(1, nRays, nLines.numLines));
                                elseif strcmp(nLines.spacing,'random')
                                    samps = randi(nRays,[nLines.numLines,1]);
                                else
                                    error('Unknown spacing parameter %s\n',nLines.spacing);
                                end
                                rays.drawSamples = samps;
                            end
                            xCoordVector = [rays.origin(samps,3) intersectPosition(samps,3) NaN([nLines.numLines 1])]';
                            yCoordVector = [rays.origin(samps,2) intersectPosition(samps,2) NaN([nLines.numLines 1])]';
                            xCoordVector = real(xCoordVector(:));
                            yCoordVector = real(yCoordVector(:));
                            figure(obj.rays.plotHandle);
                            line(xCoordVector,  yCoordVector ,'Color',lColor,'LineWidth',lWidth,'LineStyle',lStyle);
                            pause(0.1);
                        end
                    else
                        %nLines is a number 
                        % In this case, random sampling is assumed.
                        % If nLines = false, or < 0 no drawing happens
                        if (nLines > 0)
                            if (lensEl ==1)
                                rays.plotHandle = vcNewGraphWin;
                                obj.draw();
                                samps = randi(nRays,[nLines,1]);
                                rays.drawSamples = samps;
                            end
                            xCoordVector = [rays.origin(samps,3) intersectPosition(samps,3) NaN([nLines 1])]';
                            yCoordVector = [rays.origin(samps,2) intersectPosition(samps,2) NaN([nLines 1])]';
                            xCoordVector = real(xCoordVector(:));
                            yCoordVector = real(yCoordVector(:));
                            figure(rays.plotHandle);
                            pause(0.1); 
                            line(xCoordVector,  yCoordVector ,'Color',lColor,'LineWidth',lWidth,'LineStyle',lStyle);
                            
                            
                            
                            %TODO:::!!!
%                             vcNewGraphWin;
%                             for i = 1:10
%                             subplot(ceil(10/3), 3, i);
%                             plot(1:10, 1:10)
%                             end
%                             rays.plotPhaseSpace();
                        end
                       
                    end
                else
                    % This is an aperture plane.  sRadius == 0
                    intersectZ = repmat(curEl.sCenter(3), [nRays 1]); 
                    intersectT = (intersectZ - rays.origin(:, 3))./rays.direction(:, 3);
                    repIntersectT = repmat(intersectT, [1 3]);
                    intersectPosition = rays.origin + rays.direction .* repIntersectT;
                    curAperture = min(curEl.apertureD, obj.apertureMiddleD)/2;
                    
                    % Added for ppsfC aperture tracking
                    % (What does that mean?)
                    if(isa(rays, 'ppsfC'))
                        rays.aMiddleInt.XY = 0;
                        rays.aMiddleInt.XY = intersectPosition(:,1:2);  %only X-Y coords
                        rays.aMiddleInt.Z  = intersectZ;    %aperture Z
                        % passedCenterAperture = true;
                    end
                end
                
                % Set rays outside of the aperture to NaN
                outsideAperture = intersectPosition(:, 1).^2 + intersectPosition(:, 2).^2 >= curAperture^2;
                rays.origin(outsideAperture, : ) = NaN;
                rays.direction(outsideAperture, : ) = NaN;
                rays.wavelength(outsideAperture) = NaN;
                rays.waveIndex(outsideAperture) = NaN;
                intersectPosition(outsideAperture, :) = NaN;
                prevN(outsideAperture) = NaN;
                
                % Handle special case with ppsfCs
                if(isa(rays,'ppsfC'))
                    rays.aEntranceInt.XY(outsideAperture, :) = NaN;
                    rays.aMiddleInt.XY(outsideAperture, :) = NaN;                
                    rays.aExitInt.XY(outsideAperture, :) = NaN;    
                end
                
                % Apply Snell's law to the spherical surface rays.
                % Determine the new ray directions and origin
                % 
                if(curEl.sRadius ~= 0)
                    
                    %in bounds case - perform vector Snell's law
                    repCenter = repmat(curEl.sCenter, [nRays 1]);
                    normalVec = intersectPosition - repCenter;  %does the polarity of this vector matter? YES
                    normalVec = normalVec./repmat(sqrt(sum(normalVec.*normalVec, 2)),[1 3]); %normalizes each row
                    
                    %which is the correct sign convention? This is correct
                    if (curEl.sRadius < 0)  
                        normalVec = -normalVec;
                    end
                    
                    % Modify the index of refraction depending on wavelength
                    % TODO: have this be one of the input parameters (N vs. wavelength)
                    %                     if (curEl.n ~= 1)
                    %                         curN = (rays.wavelength - 550) * -.04/(300) + curEl.n;
                    %                     else
                    %                         curN = ones(length(rays.wavelength), 1);
                    %                     end
                    
                    %                     waveIndex = find([obj.wave' obj.nWave']== rays.wavelength);  %this doesn't work yet
                    %this was supposed to convert wavelength to waveIndex -
                    %but I couldn't find a way to vectorize it, so instead
                    %we precompute it
                    
                    liveIndices = ~isnan(rays.wavelength);
                    curN = ones(size(prevN));
                    curN(liveIndices) = curEl.n(rays.waveIndex(liveIndices));  %deal with nans
                    curN(~liveIndices) = nan;
                    
                    % curN = ones(length(rays.wavelength), 1) * curEl.n;
                    % Snell's law index of refraction ratios at surface
                    % boundary
                    ratio = prevN./curN;    
                    
                    % Vector form of Snell's Law
                    c = -dot(normalVec, rays.direction, 2);
                    repRatio = repmat(ratio, [1 3]);
                    
                    
                    %update the direction of the ray
                    rays.origin = intersectPosition;
     %plot phase -space for now - deal with
     %subplotting later
     
                    %plot phase space right before the lens, before the rays are bent
                    if (lensEl ==1)
                        rays.plotPhaseSpace();  %this is before the change in position
                    end
                                
                    
                    % Use bsx for speed.
                    % Simplify the line
                    newVec = repRatio .* rays.direction + repmat((ratio.*c -sqrt(1 - ratio.^2 .* (1 - c.^2))), [1 3])  .* normalVec;
                    rays.direction = normvec(newVec,'dim',2);
                    % newVec2 = newVec./repmat(sqrt(sum(newVec.*newVec, 2)), [1 3]); %normalizes each row
                    % vcNewGraphWin; plot(newVec(:),newVec2(:),'.');
                    
                    
                    
                    
                    prevN = curN;  %note: curN won't change if the aperture is the overall lens aperture
                    
                end
                
                % N.B. No need to update the direction in the case of an
                % aperture.
                
                % HURB diffraction calculation
                if (obj.diffractionEnabled)
                    obj.rtHURB(rays, intersectPosition, curEl.apertureD/2);  
                end
                
                %  xlabel('Depth (mm)');
                %  ylabel('Y field height (mm)')
                %  % iterate previous z
                % prevSurfaceZ = prevSurfaceZ + curEl.offset;
            end
        end
        
        function obj = rtIdealThroughLens(obj, rays, nLines)
            % Traces lens.rays through the ideal lens - what is that?
            % This is different from the Realistic one above, but we should
            % figure out how to make them pretty much the same.

            lWidth = 0.5; lColor = [0 0.5 1];
 
            %---------------- lens refraction code  -----
            %when intersecting ideal lens, change the direction to intersect the
            %inFocusPosition, and update the origin
            lensIntersectT = (obj.centerZ - rays.origin(:,3))./ rays.direction(:,3);
            lensIntersectPosition = rays.origin +  repmat(lensIntersectT, [1 3]) .* rays.direction;
            
            % Debug visualization
            if (nLines)
                
                samps = randi(size(rays.origin,1),[nLines,1]);
                
                xCoordVector = [rays.origin(samps,3) lensIntersectPosition(samps,3) NaN([nLines 1])]';
                yCoordVector = [rays.origin(samps,2) lensIntersectPosition(samps,2) NaN([nLines 1])]';
                xCoordVector = real(xCoordVector(:));
                yCoordVector = real(yCoordVector(:));
                line(xCoordVector,  yCoordVector ,'Color',lColor,'LineWidth',lWidth);
                pause(0.2);
                
            end
            
            %calculate new direction
            %             newRays = rayObject(); % added
            rays.origin = lensIntersectPosition;
            rays.direction = repmat(obj.inFocusPosition , [size(rays.origin,1) 1 ]) - rays.origin;
            rays.wavelength = rays.wavelength;
            
            % diffraction HURB calculation
            if (obj.diffractionEnabled)
                obj.rtHURB(rays, lensIntersectPosition, obj.apertureMiddleD/2);
            end
            %---------------- lens refraction code  -----
            
        end
    end
    
end