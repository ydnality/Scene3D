%% p_Figure1 for 2014 OSA Conference
%
%  In the ray trace calculation imagine the z-axis is positive on the right
%  and negative on the left.  The scene is located on the left.  The lens
%  is placed near zero (and has some thickness).  The film is on the right
%  (z is positive).
%
%  
% AL

%%
s_initISET

%% Initialize scene, lens and film properties

% We will loop through the lens positions
% pX = 0; pY = 0; pZ = -20000;   % millimeters

pX = -1:1; pY = -1:1; pZ = -20000:5000:-10000;   % millimeters
[X, Y, Z] = meshgrid(pX,pY,pZ);
pointSources = [X(:), Y(:), Z(:)];

% pointSources = [pX pY pZ];     % large distance test

% Create the film plane
wave = 400:100:700;            % Wavelength
wList = 1:length(wave);
fX = 0; fY = 0; fZ = 53;       % mm

% Film width and height
fW = 1;  % mm
fH = 1;  % mm

% TODO:  Get rid of wList part
% 
film = filmObject([fX fY fZ],[fW fH], wave, [wave(:) wList(:)], []);

%% Describe the lens

% At some point, we will read in a lens description file


% Multicomponent lens properties
% This goes from the light through the lens to the film
offset   = [0 1.5 1.5];   % Distances between surfaces (deltas) 
radius   = [67 0 -67];    % Radius of curvature, 0 means aperture
aperture = [4 1 4];       % Circular apertures, these are the radii in mm

% Index of refraction to the right of each surface
n = [1.67 0 1, 
     1.7 0 1,
     1.8 0 1];    

nSamples = 101;           % On the first aperture. x,y, before cropping

% May not be needed ... AL
lX = 0; lY = 0; lZ = -1.5;
lensCenterPosition = [lX lY lZ];  % Eventually calculate this given the lens file


idx = find(radius==0);  % This is the middle of the lens aperture size
fLength = 50;           % mm.  We should derive this using the lensmaker's equation
% For multiple lenses, we add up the power using something from the web

diffractionEnabled = false;
lens = lensRealisticObject(offset,radius,aperture, n(1,:), aperture(idx), fLength, lensCenterPosition, diffractionEnabled);
lens.calculateApertureSample([nSamples nSamples]);

%% Create the PSF for each point

for curInd = 1:size(pointSources, 1);
    
    %calculate the origin and direction of the rays
    rays = lens.rayTraceSourceToLens(pointSources(curInd, :));
    
    %duplicate the existing rays, and creates one for each
    %wavelength
    rays.expandWavelengths(film.wave);
    
    %lens intersection and raytrace
    lens.rayTraceThroughLens(rays);

    %intersect with "film" and add to film
    rays.recordOnFilm(film);
end



%% Loop on wavelength and depth to create PSFs 

% These psfs will be for different field heights, depths, and wavelengths


%% Show pictures

%%


%%