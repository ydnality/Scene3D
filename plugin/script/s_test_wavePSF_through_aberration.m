%s_test_wavePSF__through_aberration.m
%
%Compute the wave PSF for a psfCamera (in Scene3D). 
% The implemented method is the following:
% 1) Estimate the wavefront aberration (seidel aberration)
% 2) Exploit the Fourier Optics, for incoherent light, to compute the PSF 
%
% MP Vistasoft, 2014

%% Initialize Iset
s_initISET

%% 1.1 Declare point source in object space  and refractive indices of both object and image space

% The center of the first camera aperture is usually at [0 0 0].
% The object positions data are in -z.  We are using a right-handed
% coordinate system. 
pointSource = [10 20 -10000000];  % A very distant point.
%Refractive indices
n_ob = 1 ; %object space
n_im = 1; %image space


%% 1.2 Create a multiple lens system declaring camera properties

wave = 400:50:600;

% Declare film
filmPosition = [0 0 51.2821	];  % Good for 2Elens
% filmPosition = [0 0 37.4];        % Good for dgauss.50mm.  True focal about 37.3mm
% filmPosition = [0 0 107]; 
filmDiag = 3;  % Millimeters
filmSize = [filmDiag/sqrt(2) filmDiag/sqrt(2)];
resolution =  [300 300 length(wave)];
film = pbrtFilmC('position', filmPosition, 'size', filmSize, 'wave', wave, 'resolution', resolution);
%
lensFile = fullfile(s3dRootPath, 'data', 'lens', '2ElLens.mat');
% lensFile = fullfile(s3dRootPath, 'data', 'lens', 'dgauss.50mm.mat');
import = load(lensFile,'lens');
lens = import.lens;
lens.apertureMiddleD = 4;
lens.set('wave', wave);

% Matrix of n for each surface element.
% Apertures are 0
n = lens.get('nArray');

%% 1.3 CREATE a psfCamera structure (lightSource+lens+film [sensor])
psfCamera1 = psfCameraC('lens', lens, 'film', film, 'pointSource', pointSource);

% DO YOU WANT THE FILM IN FOCUS????  SPECIFY THE WAVELENGTH
wave0=500; %nm
psfCamera1.autofocus(wave0,'nm',n_ob,n_im);

%% 1.4 Compute and append the Black Box Model of the psfCamera 
psfCamera1.bbmCreate(n_ob,n_im);

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
[out1]=psfCamera1.drawPupil(wave0,wave,textF);
% [out1]=drawPupil(EnP,ExP,wave0,wave,textF);
%Point source
[out2]=psfCamera1.drawPoint(pointSource,wave0,wave,'y');
% Limiting rays
[out3]=psfCamera1.drawComaRay(pointSource,'entrancepupil',wave0,wave,'y');
% MARGINAL rays
% [out3a]=psfCamera1.drawMarginalRay(pointSource,'entrancepupil',wave0,wave,'y');
% and then the image point
[out4]=psfCamera1.drawPoint(imagePoint,wave0,wave,'y');
% Limiting rays
[out5]=psfCamera1.drawComaRay(imagePoint,'exitpupil',wave0,wave,'y');


%% EXAMPLE: Image a new point
pSource0=[100 100 -100000000];
[imagePoint0]=lens1.findImagePoint(pSource0,n_ob,n_im);


