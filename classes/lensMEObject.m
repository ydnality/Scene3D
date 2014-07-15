classdef lensMEObject <  handle
    % Create a multiple element lens object
    %
    %   lens = lensMEObject(elOffset, elRadius, elAperture, elN, aperture, focalLength, center);  
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
    %   lensMEObject()
    %
    % AL Vistasoft Copyright 2014
    
    properties
        name = 'default';
        type = 'multi element lens';            % We might be removing this soon (why? - BW)
        surfaceArray = lensSurfaceObject();     % Set of spherical surfaces and apertures
        diffractionEnabled = false;% Not implemented yet
        wave = 400:50:700;         % nm
        focalLength = 50;          % mm, focal length of multi-element lens
        apertureMiddleD = 1;       % mm, diameter of the middle aperture
        apertureSample = [11 11];  % Number of spatial samples in the aperture.  Use odd number
        centerZ = 0;    %theoretical center of lens (length-wise) in the z coordinate
        
        centerRay = [];   %for use for ideal Lens
        inFocusPosition = [0 0 0];  
    end
    
    
    % Private functions only called within this class
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
            % Used when we read a lens file. The function computes the center
            % position of each spherical lens element from the array of element
            % offsets and radii.
            
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
            obj.surfaceArray = lensSurfaceObject();
            
            %compute surface array centers
            centers = obj.centersCompute(sOffset, sRadius);
            
            for i = 1:length(sOffset)
                obj.surfaceArray(i) = ...
                    lensSurfaceObject('sCenter', centers(i, :), ...
                    'sRadius', sRadius(i), ...
                    'apertureD', sAperture(i), 'n', sN(i, :));
            end
            
        end
        
        function obj = rtRealisticThroughLens(obj, rays, nLines)
            % lens.rtRealisticThroughLens(rays,true/false)
            % 
            % Input:  lens object 
            %         ray structure 
            %         Number of lines for rendering
            %
            % The rays are generated by a call to lens.rtSourceToEntrance,
            % which takes a point input and generates the ray positions and
            % angles at the entrance plane.
            %
            % On return, the rays, which start out representing the rays at
            % the entrance aperture, are changed to be the position and
            % direction of the rays at the exit aperture.
            %
            % See also psfCameraObject.estimatePSF
            %
            % TODO:  Simplify this code
            %        Extend the drawing to the film plane; this should
            %        probably happen in psfCameraObject.estimatePSF.
            
            % The order is from furthest from film to film, which is also
            % how the rays pass through the optics.  We should add the film
            % draw, too.  Either here or in estimatePSF.
            
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
                if (curEl.sRadius ~= 0)
                    
                    % Spherical element
                    repCenter = repmat(curEl.sCenter, [nRays 1]);
                    repRadius = repmat(curEl.sRadius, [nRays 1]);
                    
                    % What is this formula?
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
                    
                    %added new structure to nLines to specify
                    %whether we are randomly sampling or uniformly
                    %sampling
                    
                 
                    drawLines = false;
                    
                    %struct handling - one option is to specify
                    %nLines.spacing and nLines.numLines.  .spacing has
                    %choices of 'uniform' or 'random'.  .numLines must be
                    %an integer.
                    if (isfield(nLines, 'spacing') && isfield(nLines,'numLines'))
                        if (nLines.numLines > 0)
                            if lensEl ==1
                                obj.draw();

                                if (strcmp(nLines.spacing, 'uniform'))
                                    samps = round(linspace(1, nRays, nLines.numLines));
                                else
                                    samps = randi(nRays,[nLines.numLines,1]);
                                end
                            end
                            xCoordVector = [rays.origin(samps,3) intersectPosition(samps,3) NaN([nLines.numLines 1])]';
                            yCoordVector = [rays.origin(samps,2) intersectPosition(samps,2) NaN([nLines.numLines 1])]';
                            xCoordVector = real(xCoordVector(:));
                            yCoordVector = real(yCoordVector(:));
                            line(xCoordVector,  yCoordVector ,'Color',lColor,'LineWidth',lWidth,'LineStyle',lStyle);
                            pause(0.2);
                        end
                    else
                        %number handling - another case is if nLines is a
                        %number - in this case, random sampling is assumed,
                        %and nLines = false, or 0 will assume no drawing.
                        if (nLines > 0)
                            if (lensEl ==1)
                                obj.draw();
                                samps = randi(nRays,[nLines,1]);
                            end
                            xCoordVector = [rays.origin(samps,3) intersectPosition(samps,3) NaN([nLines 1])]';
                            yCoordVector = [rays.origin(samps,2) intersectPosition(samps,2) NaN([nLines 1])]';
                            xCoordVector = real(xCoordVector(:));
                            yCoordVector = real(yCoordVector(:));
                            line(xCoordVector,  yCoordVector ,'Color',lColor,'LineWidth',lWidth,'LineStyle',lStyle);
                            pause(0.2); 
                        end
                       
                    end
                    
                    %if user has specified to draw lines, then we will.
                    if (drawLines)
                        
                    end
                    
                else
                    % This is an aperture plane.  sRadius == 0
                    intersectZ = repmat(curEl.sCenter(3), [nRays 1]); 
                    intersectT = (intersectZ - rays.origin(:, 3))./rays.direction(:, 3);
                    repIntersectT = repmat(intersectT, [1 3]);
                    intersectPosition = rays.origin + rays.direction .* repIntersectT;
                    curAperture = min(curEl.apertureD, obj.apertureMiddleD)/2;
                    
                    % Added for ppsfObject aperture tracking
                    % (What does that mean?)
                    if(isa(rays, 'ppsfObject'))
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
                
                % Handle special case with ppsfObjects
                if(isa(rays,'ppsfObject'))
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
                    newVec = repRatio .* rays.direction + repmat((ratio.*c -sqrt(1 - ratio.^2 .* (1 - c.^2))), [1 3])  .* normalVec;
                    rays.direction = normvec(newVec,'dim',2);
                    % newVec2 = newVec./repmat(sqrt(sum(newVec.*newVec, 2)), [1 3]); %normalizes each row
                    % vcNewGraphWin; plot(newVec(:),newVec2(:),'.');
                    
                    %update the direction of the ray
                    rays.origin = intersectPosition;
                    prevN = curN;  %note: curN won't change if the aperture is the overall lens aperture
                    
                end
                
                % N.B> No need to update the direction in the case of an
                % aperture.
                
                % HURB diffraction calculation
                if (obj.diffractionEnabled)
                    obj.rtHURB(rays, intersectPosition, curEl.apertureD/2);  
                end
                
                
                % iterate previous z
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
            
            %debug visualization
            if (nLines)
                
                samps = randi(size(rays.origin,1),[nLines,1]);
                
                xCoordVector = [rays.origin(samps,3) lensIntersectPosition(samps,3) NaN([nLines 1])]';
                yCoordVector = [rays.origin(samps,2) lensIntersectPosition(samps,2) NaN([nLines 1])]';
                xCoordVector = real(xCoordVector(:));
                yCoordVector = real(yCoordVector(:));
                line(xCoordVector,  yCoordVector ,'Color',lColor,'LineWidth',lWidth);
                pause(0.2);
                
                %                 for i = 1:size(rays.direction,1)
                %                     hold on;
                %                     line([rays.origin(i,3) lensIntersectPosition(i,3) ], [rays.origin(i,2) lensIntersectPosition(i,2)] ,'Color','b','LineWidth',1);
                %                 end
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
    
    methods (Access = public)
        %Multiple element lens constructor
        %TODO: error handling
        function obj = lensMEObject(varargin)
            %  surfaceList = lensReadFile(fName);
            %  lensME = lensMEObject('surfaceArray',surfList);
            
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
            %             obj.setElements(elOffset, elRadius, elAperture, elN);
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
                        obj.surfaceArray(ii).wave = val;
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
        
        
        %%
        function fileRead(obj, fullFileName)
            % Reads PBRT lens file matrix of data
            % The file has focal length information added to the header
            % This function converts the PBRT matrix of data into the format
            % that we use for setting up the multielement lens.
            
            % Open the lens file
            %fid = fopen(fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat'));
            fid = fopen(fullFileName);
            
            % Read each of the lens and close the file
            import = textscan(fid, '%s%s%s%s', 'delimiter' , '\t');
            fclose(fid);
            
            %first find the start of the lens line demarked by #   radius
            firstColumn = import{1};
            
            %Scan past the comment lines, looking for the data
            % Check the first entry in each line for #
            % When it is not, read the rest of the file as a
            % d = fread(mumble,double)
            % d= reshape(d,X,4);
            %
            % radius axpos N aperture
            % The current test is just to see if the word radius is on the
            % line.  This might be improved.
            continu = true;
            i = 1;
            while(continu && i <= length(firstColumn))
                compare = regexp(firstColumn(i), 'radius');
                if(~(isempty(compare{1})))
                    continu = false;
                end
                i = i+1;
            end
            % i is the row where the data begin
            
            % The next lens are the matrix data
            % put data into lens object
            radius = str2double(import{1});
            radius = radius(i:length(firstColumn));
            
            
            % Change from pbrt Scene3D format to raytrace Scene3D format
            % In PBRT, the row has the offset from the previous surface.  In
            % PBRT the data are read from the bottom up.  The last row has no
            % offset.
            % In PBRT, we trace from the sensor to the scene.
            % In Scene3D we trace from the scene to the sensor.
            % So, the offsets are shifted down.  This means:
            %
            offset = str2double(import{2});
            offset = offset(i:length(firstColumn));
            offset = [0; offset(1:(end-1))];
            
            % Index of refraction in the 3rd column
            N = str2double(import{3});
            N = N(i:length(firstColumn));
            
            % Diameter of the aperture (or maybe radius.  Must determine).
            aperture = str2double(import{4});
            aperture = aperture(i:length(firstColumn));
            
            %modify the object and reinitialize
            obj.elementsSet(offset, radius, aperture, N);
        end
        
        %%
        function numEls = numEls(obj)
            % This is now lens.get('n surface')
            warning('Use lens.get(''n surface'')');
            numEls = length(obj.surfaceArray);
        end
        
        %% Draw the spherical and aperture components
        function obj =  draw(obj)
            % Draw the the multi-element lens in a new graph window. 
            % Helpful for debugging
            
            % Create the figure and set the parameters
            vcNewGraphWin([],'wide');
            axis equal;
            lWidth = 2; lColor = 'k';  % Drawing parameters

            % We draw one surface/aperture at a time
            nSurfaces = obj.get('n surfaces');
            for lensEl = 1:nSurfaces
                
                % Get the current surface element
                curEl = obj.surfaceArray(lensEl);
                
                if (curEl.sRadius ~= 0)
                    %% Simpler code
                    testme = false;
                    if testme
                        [z,y] = circlePoints([curEl.sCenter,0],curEl.sRadius);
                        if curEl.sRadius > 0
                            l = (abs(y) < curEl.apertureD/2) & (z < curEl.sCenter(3));
                        else
                            l = (abs(y) < curEl.apertureD/2) & (z > curEl.sCenter(3));
                        end
                        p = plot(z(l),y(l),'k.'); set(p,'markersize',1); hold on;
                    else
                        % Draw a spherical element
                        
                        % Get previous and next elements, within bounds
                        nextEl = obj.surfaceArray(min(lensEl+1, end));
                        prevEl = obj.surfaceArray(max(lensEl-1, 1));
                        
                        % Lens elements do NOT always end when the neighboring
                        % element begins.  this allows for a fudge factor.  This
                        % won't matter too much because the aperture radius will
                        % be the limiting factor. (I don't understand this
                        % comment, other than this is a fudge factor somehwere.
                        % BW).
                        delta = 10;
                        
                        % Determine the minimum/maximum z-positions for the
                        % curve. Which side has the position depends on the
                        % sign of the curvature
                        if (curEl.sRadius > 0 )
                            % The fudge factor is used here
                            leftBoundary  = curEl.get('zIntercept');
                            rightBoundary = nextEl.get('zIntercept') + delta;
                        else
                            % Here is fudge again
                            leftBoundary  = prevEl.get('zIntercept') - delta;
                            rightBoundary = curEl.get('zIntercept');
                        end
                        
                        % This is the range of z values we will consider.
                        zPlot = linspace(leftBoundary, rightBoundary, 100);
                        
                        
                        % Solve for the points on the curve for this lens
                        % surface element.  We are drawing in the z-y plane
                        % because the z-axis is the horizontal axis, and the
                        % y-axis is the vertical. The center of the sphere is
                        % at (0,0,z), so the formula is
                        %
                        %     r^2 = (x)^2 + (y)^2 + (z - c)^2
                        %
                        % But since we are in the x=0 plane this simplifies to
                        %
                        %   r^2 = (y)^2 + (z - c)^2
                        %
                        % We solve for y in terms of z.
                        % We get the positive and negative y-values
                        yPlot  =  sqrt(curEl.sRadius^2 - (zPlot - curEl.sCenter(3)) .^2);
                        yPlotN = -sqrt(curEl.sRadius^2 - (zPlot - curEl.sCenter(3)) .^2);
                        
                        % NOTE: We may have problems with a concave lens
                        % because of how we are choosing the z range - but
                        % these are rare.
                        %
                        % We will plot for the range of y values that are less
                        % than the radius of the spherical surface.
                        withinRange = (yPlot < curEl.apertureD/2);
                        
                        % The positive solutions
                        l = line(zPlot(withinRange), yPlot(withinRange));
                        set(l,'linewidth',lWidth,'color',lColor);
                        
                        % The negative solutions
                        l = line(zPlot(withinRange), yPlotN(withinRange));
                        set(l,'linewidth',lWidth,'color',lColor);
                    end
                    
                else
                    %Draw the aperture opening if radius = 0
                    
                    %TODO: draw the difference between specified aperture
                    %from file and specified aperture from object instance
                    
                    %right now: take the minimum value
                    curAperture = min(curEl.apertureD/2, obj.apertureMiddleD/2);
                    
                    l = line(curEl.sCenter(3) * ones(2,1), -1*[curEl.apertureD/2 curAperture]);
                    set(l,'linewidth',lWidth,'color',lColor);
                    l = line(curEl.sCenter(3) * ones(2,1), [curAperture curEl.apertureD/2]);
                    set(l,'linewidth',lWidth,'color',lColor);

                end

            end
        end
        

        
        
        function obj = rtThroughLens(obj, rays, nLines, rtType)
            % lens.rtThroughLens(rays,true/false)
            % 
            % Inputs:  lens object,rays at the entrance aperture,
            % calculates the rays at the exit aperture.  The exiting rays
            % are represented by the 
            %
            % On return, the input variable rays, which starts out
            % representing the rays at the entrance aperture, is changed to
            % be the position and direction of the rays at the exit
            % aperture.
            %
            % Simplify the code.
            % Also, why is there similar code in ppsfCamera?
            
            % The order is from furthest from film to film, which is also
            % how the rays pass through the optics.
            
            if (ieNotDefined('nLines')), nLines = false; end
            if (ieNotDefined('rtType')), rtType = 'realistic'; end
             rtType = ieParamFormat(rtType);
             switch rtType
                 case 'ideal'
                     obj = obj.rtIdealThroughLens(rays, nLines);
                 case 'realistic'
                     obj = obj.rtRealisticThroughLens( rays, nLines);
                 case 'linear'
                     error ('not implemented yet');
                 otherwise
                     error ('unknown ray trace type');
             end
        end
            

        function rays = rtSourceToEntrance(obj, pointSource, ppsfObjectFlag, jitterFlag, rtType)
            % Ray trace from a point to the aperture grid on the first lens
            % surface (the one furthest from the sensor).  
            %  
            % Rays have 1 wavelength assigned to them.  There is no
            % particular wavelength dependence in air, so there is no need
            % to have multiple indices of refraction or wavelength at this
            % moment.
            %
            % The rays will be "expanded" out later to save on computations
            % to handle wavelength differences.
            %
            % The object is a lensMEObject
            %     pointSource    -  a 3 vector
            %     ppsfObjectFlag -  specifies whether we are computing the
            %     plenoptic PSF (pPSF) or not.  Default: false (PSF only)
            %
            % See also: rtEntranceToExit 
            
            if (ieNotDefined('ppsfObjectFlag')), ppsfObjectFlag = false;end
            if (ieNotDefined('jitterFlag')), jitterFlag = false; end
            if (ieNotDefined('rtType')), rtType = 'realistic'; end
            rtType = ieParamFormat(rtType);
            
            % Define rays object
            if (~ppsfObjectFlag), rays = rayObject();
            else                  rays = ppsfObject();
            end
            
            
            %deal with different ray trace types
            switch rtType
                case 'realistic'
                    frontRTSurface = -obj.get('totalOffset');
                case 'ideal'
                    frontRTSurface = obj.centerZ;
                    
                    % --- center ray calculation-------------
                    
                    %trace ray from point source to lens center, to image.  This helps
                    %determine the point of focus
                    obj.centerRay.origin = pointSource;
                    
                    % This could be a call to rayDirection
                    obj.centerRay.direction = [ 0 0 obj.centerZ] - obj.centerRay.origin;
                    obj.centerRay.direction = obj.centerRay.direction./norm(obj.centerRay.direction);
                    
                    %calculate the z-position of the in-focus plane using
                    %thin lens equation. 
                    inFocusDistance = 1/(1/obj.focalLength - -1/pointSource(3));
                    
                    % Calculates the 3-vector for the in-focus position.
                    % The in-focus position is the intersection of the
                    % in-focus plane and the center-ray
                    inFocusT = (inFocusDistance - obj.centerRay.origin(3))/obj.centerRay.direction(3);
                    obj.inFocusPosition = obj.centerRay.origin + inFocusT .* obj.centerRay.direction;
                    % --------center ray calculation -------
                case 'linear'
                    error('not implemented yet')
                otherwise
                    error('unknown ray trace method:  %s\n',rtType)
            end
             
            
            % Create rays from the point source to each aperture grid point
            % The new origin is the position of the current point source
            % The jitterFlag scatters the points on the front surface a bit
            % to improve the rendering.
            aGrid   = obj.apertureGrid(jitterFlag);
            
            % Set Z to the position of the front surface
            nPts    = numel(aGrid.X(:));
            aGrid.Z = repmat(frontRTSurface,[nPts,1]);
            
            % These are the end points of the ray in the aperture plane
            ePoints = [aGrid.X(:),aGrid.Y(:),aGrid.Z(:)];  
            
            % These are the directions from the point source to the end
            % points in the aperture
            rays.origin    = repmat(pointSource, [nPts, 1, 1] );
            rays.direction = rayDirection(rays.origin,ePoints);
            
            if(isa(rays,'ppsfObject'))
                % If the rays are a plenoptic point spread, AL thinks we
                % should store the XY positions at the front aperture,
                % middle aperture, and exist aperture.  The slots for this
                % information are created and initialized here.
                %
                % BW isn't sure why this is useful, but he believes it.
                rays.aEntranceInt.XY = zeros(length(aGrid.X), 2);
                rays.aEntranceInt.XY(:,1) = aGrid.X(:);
                rays.aEntranceInt.XY(:,2) = aGrid.Y(:);
                rays.aEntranceInt.Z = aGrid.Z(1);
                
                % The intersection positions (XY) of the rays in the middle
                % aperture and the exit aperture will be stored here.
                rays.aMiddleInt.XY = zeros(length(aGrid.X), 2);
                rays.aExitInt.XY = zeros(length(aGrid.X), 2);
                
            end
            
        end
        
        function obj = rtHURB(obj, rays, lensIntersectPosition, curApertureRadius)
            %Performs the Heisenburg Uncertainty Ray Bending method on the
            %rays, given a circular aperture radius, and lens intersection
            %position This function accepts both vector forms of inputs, or
            %individual inputs
            %
            % Look for cases when you can use: bsxfun ...
            %
            % This code is not readable yet.  Let's figure out the steps
            % and write clarifying functions.
            %
            
            
            ipLength = sqrt(sum(dot(lensIntersectPosition(:, (1:2)), lensIntersectPosition(:, (1:2)), 2), 2));
            
            %calculate directionS which is ....
            directionS = [lensIntersectPosition(:, 1) lensIntersectPosition(:,2) zeros(length(lensIntersectPosition), 1)];

            % And the orthogonal directionL
            directionL = [-lensIntersectPosition(:,2) lensIntersectPosition(:,1) zeros(length(lensIntersectPosition), 1)];
            
            normS = repmat(sqrt(sum(dot(directionS, directionS, 2), 2)), [1 3]);
            normL = repmat(sqrt(sum(dot(directionL, directionL, 2), 2)), [1 3]);
            divideByZero = sum(normS,2) ==0;
            
            directionS(~divideByZero, :) = directionS(~divideByZero, :)./normS(~divideByZero, :);
            directionL(~divideByZero, :) = directionL(~divideByZero, :)./normL(~divideByZero, :);
            directionS(divideByZero, :) = [ones(sum(divideByZero == 1), 1) zeros(sum(divideByZero == 1), 2)];
            directionL(divideByZero, :) = [zeros(sum(divideByZero == 1), 1) ones(sum(divideByZero == 1), 1) zeros(sum(divideByZero == 1), 1)];
            
            %                 if (norm(directionS)~= 0)
            %                     directionS = directionS./norm(directionS);
            %                     directionL = directionL./norm(directionL);
            %                 else
            %                     %this is the case that the ray hits the absolute center
            %                     %Then, there is a special case to avoid division by 0
            %                     directionS = [1 0 0];
            %                     directionL = [0 1 0];
            %                 end
            
            pointToEdgeS = curApertureRadius - ipLength;   %this is 'a' from paper  //pointToEdgeS stands for point to edge short
            pointToEdgeL = sqrt((curApertureRadius* curApertureRadius) - ipLength .* ipLength);  %pointToEdgeS stands for point to edge long
            
            lambda = rays.wavelength * 1e-9;  %this converts lambda to meters
            sigmaS = atan(1./(2 * pointToEdgeS *.001 * 2 * pi./lambda));  %the .001 converts mm to m
            sigmaL = atan(1./(2 * pointToEdgeL * .001 * 2 * pi./lambda));
            
            %this function regenerates a 2D gaussian sample and
            %returns it randOut
            %gsl_ran_bivariate_gaussian (r, sigmaS, sigmaL, 0, noiseSPointer, noiseLPointer);    %experiment for now
            [randOut] = randn(length(sigmaS),2) .* [sigmaS sigmaL];
            
            %calculate component of these vectors based on 2 random degrees
            %assign noise in the s and l directions according to data at these pointers
            noiseS = randOut(:,1);
            noiseL = randOut(:,2);
            
            %project the original ray (in world coordinates) onto a new set of basis vectors in the s and l directions
            projS = (rays.direction(: , 1) .* directionS(: ,1) + rays.direction(: , 2) .* directionS(:,2))./sqrt(directionS(:,1) .* directionS(:,1) + directionS(:,2) .* directionS(:,2));
            projL = (rays.direction(: , 1) .* directionL(:, 1) + rays.direction(: , 2 ) .* directionL(:,2))./sqrt(directionL(:,1) .* directionL(:,1) + directionL(:,2) .* directionL(:,2));
            thetaA = atan(projS./rays.direction(: , 3));   %azimuth - this corresponds to sigmaS
            thetaE = atan(projL./sqrt(projS.*projS + rays.direction(:, 3).* rays.direction(: , 3)));   %elevation - this corresponds to sigmaL
            
            %add uncertainty
            thetaA = thetaA + noiseS;
            thetaE = thetaE + noiseL;
            
            %convert angles back into cartesian coordinates, but in s,l space
            newprojL = sin(thetaE);
            smallH = cos(thetaE);   %smallH corresponds to the projection of the ray onto the s-z plane
            newprojS = smallH .* sin(thetaA);
            rays.direction(:, 3) = smallH .* cos(thetaA);
            
            %convert from s-l space back to x-y space
            rays.direction(:, 1) = (directionS(:, 1) .* newprojS + directionL(:, 1) .* newprojL)./sqrt(directionS(:,1) .* directionS(:,1) + directionL(:,1) .* directionL(:,1));
            rays.direction(:, 2) = (directionS(:, 2) .* newprojS + directionL(:, 2) .* newprojL)./sqrt(directionS(:,2) .* directionS(:,2) + directionL(:,2) .* directionL(:,2));
            normDirection = repmat(sqrt(sum(dot(rays.direction, rays.direction, 2),2)), [1 3]);
            rays.direction = rays.direction./normDirection;
            
            %reassign ray
            %                 rays.origin(i,:) = curRay.origin;
            %                 rays.direction(i, :) = curRay.direction;
            %                 rays.wavelength(i,:) = curRay.wavelength;
            %             end
        end
        
        
    end
    
end