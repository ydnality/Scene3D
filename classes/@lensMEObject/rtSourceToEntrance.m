function rays = rtSourceToEntrance(obj, pointSource, ppsfObjectFlag, jitterFlag, rtType)
% Ray trace from a point to the aperture grid on the first lens surface
% (the one furthest from the sensor).
%
%  rays = rtSourceToEntrance(obj, pointSource, ppsfObjectFlag, jitterFlag, rtType)
%
% Rays have 1 wavelength assigned to them.  There is no particular
% wavelength dependence in air, so there is no need to have multiple
% indices of refraction or wavelength at this moment.
%
% The rays will be "expanded" later to save on computations to
% handle wavelength differences.
%
% The object is a lensMEObject
%     pointSource    -  a 3 vector
%     ppsfObjectFlag -  specifies whether we are computing the
%     plenoptic PSF (pPSF) or not.  Default: false (PSF only)
%
% See also: rtEntranceToExit

if ieNotDefined('ppsfObjectFlag'), ppsfObjectFlag = false; end
if ieNotDefined('jitterFlag'),     jitterFlag = false;     end
if ieNotDefined('rtType'),         rtType = 'realistic';   end

% Define rays object
if (~ppsfObjectFlag), rays = rayObject();
else                  rays = ppsfObject();
end

%deal with different ray trace types
rtType = ieParamFormat(rtType);
switch rtType
    case 'realistic'
        frontRTSurface = -obj.get('totalOffset');
    case 'ideal'
        frontRTSurface = obj.centerZ;
        
        % --- center ray calculation-------------
        
        %trace ray from point source to lens center, to image.
        %This helps determine the point of focus
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
        error('Linear method not implemented yet')
    otherwise
        error('Unknown ray trace method:  %s\n',rtType)
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