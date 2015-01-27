%% s_s3dTransverseCATemp
%
% Show how changing the aperture position with respect to the lens causes a
% shift in the magnification (transverse chromatic aberration) with respect
% to wavelength
%
% This is built on the notes DHB left us about 18 months ago
%
% We should clean this up for demo'ing and check it with the other one that
% MP is looking over.
%
% AL/BW Vistasoft Team, Copyright 2014

%%
s_initISET

%% Make a lens with an aperture in the middle.
lens = lensC;
lens.draw;

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

film = filmC;
% Make a point source
ps = [-5 0 -100];

ppsfCamera = ppsfCameraC('lens',lens,'film',film,'point source',ps);

%% 
nLines = 0; jitterFlag = true;

% Should we add to the psf or should we start fresh?  We need to be
% clearer.
ppsfCamera.estimatePSF(nLines,jitterFlag);

oi = ppsfCamera.oiCreate;
vcAddObject(oi); oiWindow;

%% Shifting the aperture position to create transverse CA

s = lens.get('surface array',2);
p = s.get('zpos');
s.set('zpos',p - 22)
lens.sortSurfaceOrder;

nLines = 0; jitterFlag = true;

% Should we add to the psf or should we start fresh?  We need to be
% clearer.
ppsfCamera.estimatePSF(nLines,jitterFlag);

oi = ppsfCamera.oiCreate;
vcAddObject(oi); oiWindow;


%% End