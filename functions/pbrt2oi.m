function oi = pbrt2oi(fname,inputPbrt)
%Convert pbrt multispectral irradiance data into an ISET optical image
%
%   oi = pbrt2oi(fname,inputPbrt)
%
% This routine is relied upon by the main rendering functions s3dRenderOI
% and so forth.
%
% fname is a file produced by pbrt.
% inputPbrt is a pbrt object that contains information about the film and
% lens
% 
% See also:  s3dRenderOIAndDepthMap, s3dRenderOI
%
% (c) Stanford VISTA Team 2012

if ieNotDefined('fname'), error('File name required.'); end
if ieNotDefined('inputPbrt'), inputPbrt = []; end

%open file
fID = fopen(fname,'r','l');

%load size information
[A, cnt] = fscanf(fID,'%d %d %d\n',[3 1]); %#ok<NASGU>

% The FOV and lens information now comes from elsewhere.
%load lens and field of view and information
[FOV, cnt2] = fscanf(fID,'%f %f %f\n',[3 1]); %#ok<ASGLU,NASGU>

% We don't think these numbers are right yet.
% Making them right could involve PBRT hacking, or using the lens
% information that is part of the pbrtObject in Matlab.
% if (~isempty(FOV))
%     % focalLength = FOV(1)*1e-3;   %do something with this information in the future
%     % aperture    = FOV(2)*1e-3;
%     % fieldOfView = FOV(3);
% else
%     disp('no lens information!!')
% end

%Load the stored photons produced by AL's pbrt code
photons = fread(fID,prod(A),'double');

%A(2) = rows, A(1) = columns
photons = reshape(photons,A(2),A(1),A(3));
fclose(fID);

% Set the OI data
oi = oiCreate;
oi = initDefaultSpectrum(oi);

% Put the irradiance in
% Check about this 32nd, and wavelength and all that!
% For now, we always run at 400:10:700 really, but AL needed to add one
% more for something about the PBRT calculation.  Here we toss the last
% wavelength.
oi = oiSet(oi,'photons',single(photons(:,:,1:31)));

% Set the optics parameters from somewhere
% oi = oiSet(oi,'optics focal length',focalLength);
% oi = oiSet(oi,'optics fnumber',focalLength/aperture);
% oi = oiSet(oi,'fov',fieldOfView);

oi = oiSet(oi, 'photons', oiGet(oi,'photons') * 10^13);  %some normalization issues


% Calculate field of view
if isequal(class(inputPbrt),'pbrtObject')
    
    lens = inputPbrt.camera.lens;
    
    if isa(lens,'pbrtLensRealisticObject')
        oi = oiSet(oi,'optics name',lens.specFile);
        oi = oiSet(oi,'optics fnumber',lens.filmDistance/lens.apertureDiameter);       
    elseif isa(lens,'pbrtLensPinholeObject')
        % Pinholes have no real aperture size.  So, we set the f-number
        % really big.
        oi = oiSet(oi,'optics name','pinhole');
        oi = oiSet(oi,'optics fnumber',999);
    elseif isa(lens,'pbrtLensIdealObject')
        % This case is a diffraction limited lens but with an aperture of a
        % real size.
        oi = oiSet(oi,'optics name','diffraction limited');
        oi = oiSet(oi,'optics fnumber',lens.filmDistance/lens.apertureDiameter);
    end
    
    oi = oiSet(oi,'optics focal length',lens.filmDistance*1e-3);

    % Compute the horizontal field of view
    fdiag = lens.filmDiag;
    dist  = lens.filmDistance;
    x     = inputPbrt.camera.film.xresolution;
    y     = inputPbrt.camera.film.yresolution;
    d     = sqrt(x^2+y^2);  % Number of samples along the diagonal
    fwidth= (fdiag/d)*x;    % Diagonal size by d gives us mm per step
                            % multiplying by x gives us the horizontal mm  
    % Calculate angle in degrees
    fov = 2*atan2d(fwidth/2,dist);
    
elseif ~isempty(inputPbrt) && isnumeric(inputPbrt)
    % This should never happen!  Alert the user that we are completely
    % making this up.
    dist = inputPbrt;  % Treat the focal length as the location of the sensor
    fwidth = 30;      %  Thirty five millimeter diagonal, 
    fov = 2*atan2d(fwidth/2,dist);
    fprintf('inputPBRT is numeric. Arbitrarily setting the film width to %.2f mm',fwidth);
else
    % Should never get here, please.  But we do now.
    fov = 10;
    fprintf('Fov arbitrarily set to 10 deg.\n');
end

% Store the horizontal field of view in degrees in the oi
oi = oiSet(oi,'fov', fov);

return
