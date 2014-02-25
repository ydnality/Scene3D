%% Ray-tracing for an ideal lens
%
% An ideal thin lens is defined as one that obey's the thin lens equation
% (1/s1 + 1/s2 = 1/f).  Therefore, Snell's Law is not used in this script
% at all.  Instead, the direction that the rays bend are determined by the
% thin lens equation.  At each point source, a "center ray" is shot at the
% center of the lens.  The intersection of this ray and the focal-plane as
% defined by the thin lens equation determines the point of focus of the
% lens.  All other rays that are shot at the edge of the aperture will then
% intersect this ideal focal point.  
% All units are in mm.
%
%  TODO:
%   We can't use sensor here and sensorCreate in same way.
%   Somehow we should integrate this sensor stuff with ISET.
%   I don't see the lens part here.  Need more comments about what this is
%   doing, and what it tests.
%
% See also: s_3dRayTrace*.m
%
% AL Vistalab 2014


%% point sources
% declare point sources in world space.  The camera is usually at [0 0 0],
% and pointing towards -z.  We are using a right-handed coordinate system.

% [XGrid YGrid] = meshgrid(-4000:1000:4000,-4000:1000:4000);
% pointSources = [XGrid(:) YGrid(:) ones(size(XGrid(:))) * -20000];
pointSources = [0 0 -20000];


%% camera and film properties 
diffractionEnabled = true;
film = filmObject([0 0 50],[1 1], 400:10:700, [(400:10:700)' (1:31)'], []);   %(position, size,  wave, waveConversion, resolution)
lens = lensIdealObject(.1, 50, [0 0 0], diffractionEnabled);
lens.calculateApertureSample([20 20]);

vcNewGraphWin; 

%% loop through all point sources
for curInd = 1:size(pointSources, 1);
    rays = rayObject;
    %calculate the origin and direction of the rays
    lens.rayTraceSourceToLens(pointSources(curInd, :), rays);
    
    %duplicate the existing rays, and creates one for each
    %wavelength
    rays.expandWavelengths(film.wave);

    %lens intersection and raytrace
    lens.rayTraceThroughLens(rays, pointSources(curInd, :));

    % intersect with "film" and add to film
    rays.recordOnFilm(film);
end

%%
%assign as optical image
oi = oiCreate;
oi = initDefaultSpectrum(oi);
oi = oiSet(oi,'photons',film.image);
vcAddAndSelectObject(oi); oiWindow;

%%
