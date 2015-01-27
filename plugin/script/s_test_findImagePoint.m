% s_test_findImagePoint.m
%
% Test the conversion of a system from SCENE 3D to format of PSF3D
%
% MP Vistasoft, 2014

%% Initialize Iset
s_initISET

%% 1.1 Declare point source in object space

% The center of the first camera aperture is usually at [0 0 0].
% The object positions data are in -z.  We are using a right-handed
% coordinate system. 
p = psCreate(0,0,-50);
pointSource = p{1};  % A very distant point.

%% Declare film and lens 

% lensFile = fullfile(s3dRootPath, 'data', 'lens', '2ElLens.mat');
lensFile = fullfile(s3dRootPath, 'data', 'lens', 'dgauss.50mm.mat');

% Loading this twice works, but the first time it fails.  No idea why. (BW)
load(lensFile,'lens');
lens.apertureMiddleD = 4;
wave = lens.get('wave');

% Declare film.  Match aspects of it to the lens
% filmPosition = [0 0 51.2821	];  % Good for 2Elens
filmPosition = [0 0 37.4];        % Good for dgauss.50mm.  True focal about 37.3mm
filmDiag = 3;  % Millimeters
filmSize = [filmDiag/sqrt(2) filmDiag/sqrt(2)];
resolution =  [300 300 length(wave)];
film = filmC('position', filmPosition, 'size', filmSize, 'wave', wave, 'resolution', resolution);

% Matrix of n for each surface element.
% Apertures are 0
% n = lens.get('nArray');

%% Create an imaging system with a black box model

psfCamera1 = psfCameraC('lens', lens, 'film', film, 'pointSource', pointSource);
psfCamera1.bbmCreate;
disp(psfCamera1)
disp(psfCamera1.BBoxModel)

%
n_ob = 1; n_im = 1.3;  % Like air to eye interface
psfCamera1.bbmCreate(n_ob,n_im);
disp(psfCamera1)
disp(psfCamera1.BBoxModel)


%% Add the Black Box Model to the objects (lens or psfCamera or 'all')

% These three functions convert the Scene3D objects designed around the
% light field geometric representation into more complex objects that have
% wavefront information.

% This one takes a lens with a light field description and fills in the
% 'black box model (BBoxModel)' to the lens structure.  The black box model
% is a simplification of the multi-element lens.
%   
% There is a bbmGetValue (but we should probably just use 
% lens.get('bbm XXX') and lens.set('bbm XXX)

%Refractive indices
n_ob = 1 ; %object space
n_im = 1; %image space

% Add the BBoxModel to the lens
% [lens1]      = paraxAnalyzeScene3DSystem('lens',lens);
lens1 = lens;
lens1.bbmCreate(n_ob,n_im);
disp(lens1)

% ALTERNATIVE 
bbm_lens1 = lens1.bbmCreate(n_ob,n_im);
disp(bbm_lens1)


% GET OPTICAL SYSTEM STRUCTURE (used at data source for Black Box Model Field)
OptSyst1 = lens1.bbmCreate(n_ob,n_im);

% Add the BBoxModel to the camera
% psfCamera1.bbmCreate(n_ob,n_im);
% ALTERNATIVE 
bbm_psfC1 = psfCamera1.bbmCreate(n_ob,n_im);
% GET OPTICAL SYSTEM STRUCTURE ( used at data source for Black Box Model Field)
% ImgSyst1 = psfCamera1.bbmCreate(n_ob,n_im);

% Create a camera throught the BBoxModel of lens
% That's weird (BW).
% [psfCamera2]    = paraxAnalyzeScene3DSystem('all',lens,film,pointSource);
psfCamera2 = lens1.bbmCreate('all',pointSource,film,n_ob,n_im);

%% GET SEVERAL FIELDs
%
% H is for principal
% N is for a nodal point
% F is for a focal point
%
% The prime character (') is usually appended for conjugate points
%

Hobj = psfCamera1.get('bbm','object principal point');       %principal point in the object space
Him  = psfCamera1.get('bbm','image principal point');        %principal point in the image space
focalLength = psfCamera1.get('bbm','effective focal length');%effective focal length
magn        = psfCamera1.get('bbm','lateral magnification'); %lateral magnification
abcd  = psfCamera1.get('bbm','abcd matrix');                 %abcd matrix coeffs

%% FIND the image point

%
% These are lambda x position, so each row is one wavelength
w = psfCamera1.get('wavelength');
p = psfCamera1.get('bbm','gaussian image point');
[w,p]

% This is the location of the in focus image point given the point source
% position and the lens object.  You can use this value to place the film
% for perfect focus for this point.
imagePoint = lens1.findImagePoint(pointSource,n_ob,n_im);

% You can also find the in focus image point from the camera object this
% way.  This call uses the same function.
imagePoint2 = psfCamera1.get('bbm','gaussian image point')

imagePoint = imagePoint2

%% PLOT the Entrance and Exit Pupil
%
% EnP=psfCamera1.EntrancePupil;
% ExP=psfCamera1.ExitPupil;
% 
% overDiam=[EnP.diam,ExP.diam];
% maxDim=2*max(max(overDiam)); %twice the max pupil diameter
% 
wave0=550; %nm   select a wavelengt
% 
% % set text box
textF='true';

% lens.draw
psfCamera1.draw

% Graph pupil
% [out1]=psfCamera1.drawPupil(wave0,wave,textF);
% %Point source
[out2]=psfCamera1.drawPoint(pointSource,wave0,wave,'y','Point Source');
% % Limiting rays
[out3]=psfCamera1.drawComaRay(pointSource,'entrancepupil',wave0,wave,'y');
% % MARGINAL rays
% % [out3a]=psfCamera1.drawMarginalRay(pointSource,'entrancepupil',wave0,wave,'y');
% % and then the image point
[out4]=psfCamera1.drawPoint(imagePoint,wave0,wave,'y','Image Point');
% % Limiting rays
[out5]=psfCamera1.drawComaRay(imagePoint,'exitpupil',wave0,wave,'y');

% END

