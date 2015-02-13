%% s_s3dTransverseCATemp
%
% Show how changing the aperture position with respect to the lens causes a
% shift in the magnification with respect to wavelength (transverse
% chromatic aberration). 
% This is built on the notes DHB left us about 18 months ago
%
% We should clean this up for demo'ing and check it with the other one that
% MP is looking over.
%
% AL/BW Vistasoft Team, Copyright 2014

%%
ieInit

%% Make a lens with an aperture in the middle.
lens = lensC;
% lens.draw;

%% Set index of refraction

% Set the n values for all the lenses to range from 1.3 to 1.7
% Should the aperture (n = 6) have all ones?
% This breaks the MP's code.  Commenting it out lets the rest of the code
% run coorectly.
nSurfaces = lens.get('n surfaces');
for ii=1:(nSurfaces-1)
    if lens.surfaceArray(ii).sRadius ~= 0
        lens.surfaceArray(ii).n = linspace(1.3,1.5,lens.get('nwave'));
    end
end

%% Make the sensor surface

film = filmC;

%% Make a point source
ps = [-5 0 -100];

ppsfCamera = ppsfCameraC('lens',lens,'film',film,'point source',ps);

%% 

% Should we add to the psf or should we start fresh?  We need to be
% clearer.
nLines = 0; jitterFlag = true;
ppsfCamera.estimatePSF(nLines,jitterFlag);

oi = ppsfCamera.oiCreate;
vcAddObject(oi); oiWindow;

%% Shift the aperture position to create transverse CA

% The second element of the surface array is the aperture (diaphragm)
aSurface = lens.get('aperture');  % The surface with the aperture
s        = lens.get('surface array',aSurface);
pOrig    = s.get('zpos');     % Original position of the aperture

%% Show the changing magnification with aperture z position
d =[0, -30];

for ii=1:length(d)
    aSurface = lens.get('aperture');  % The surface with the aperture
    s        = lens.get('surface array',aSurface);

    s.set('zpos',pOrig + d(ii));  % Change the z position of the aperture
    lens.sortSurfaceOrder;        % Make sure the s array is updated with the new position
    % vcNewGraphWin; lens.draw;
    
    % Should we add to the psf or should we start fresh?  We need to be
    % clearer.
    nLines = 0; % No debug lines
    jitterFlag = true;
    ppsfCamera.estimatePSF(nLines,jitterFlag);
    
    p = s.get('z pos'); 
    
    % Show it in ISET
    oi = ppsfCamera.oiCreate;
    oi = oiSet(oi,'name',sprintf('Pos = %.1f',p));
    vcAddObject(oi); oiWindow;
    ppsfCamera.film.clear();
    % Plots
    %     x = 1; y = 100; w = 7;
    %     uData = plotOI(oi,'irradiance hline',[x,y]);  % (x,y) position
    %     irrad(:,ii) = uData.data(w,:)';
end
% plot(irrad)


%% Using bbm 
% ppsfCamera.bbmCreate;
% 
% ppsfCamera.bbmGetValue('effective focal length')


%% End