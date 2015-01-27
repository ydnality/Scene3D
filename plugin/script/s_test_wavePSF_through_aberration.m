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

% Specify the point source

% Angulare eccentricity (in degree)
angle_field=5; %�
% Distance 
ps_dist=0.25*1e4; %distance

% Create point source coordinate
ps_x=ps_dist.*tan(angle_field/180*pi);
pointSource = [ps_x 0 -ps_dist];  % A very distant point.

%Refractive indices
n_ob = 1 ; %object space
n_im = 1; %image space


%% 1.2 Create a multiple lens system declaring camera properties

wave = 400:50:600;

% Declare film
% filmPosition = [0 0 51.2821	];  % Good for 2Elens
filmPosition = [0 0 37.4];        % Good for dgauss.50mm.  True focal about 37.3mm
% filmPosition = [0 0 107]; 
filmDiag = 3;  % Millimeters
filmSize = [filmDiag/sqrt(2) filmDiag/sqrt(2)];
resolution =  [300 300 length(wave)];
film = filmC('position', filmPosition, 'size', filmSize, 'wave', wave, 'resolution', resolution);
%
% lensFile = fullfile(s3dRootPath, 'data', 'lens', '2ElLens.mat');
lensFile = fullfile(s3dRootPath, 'data', 'lens', 'dgauss.50mm.mat');
import = load(lensFile,'lens');
lens = import.lens;
lens.apertureMiddleD = 4;
lens.set('wave', wave);

% Matrix of n for each surface element.
% Apertures are 0
n = lens.get('nArray');

%% 1.3 CREATE a psfCamera structure (lightSource+lens+film [sensor])
psfCamera1 = psfCameraC('lens', lens, 'film', film, 'pointSource', pointSource);


bbmOLD=psfCamera1.bbmCreate(n_ob,n_im);
defocusOLD=psfCamera1.get('bbm','defocus','#wave');
paOLD=psfCamera1.get('bbm','primaryaberration','#wave');


%% GET the PSF
nSample=256; % power of 2 for Fast Fourier Transform
ntime=16;  % window width= ntime*normalized_ExPdiameter
[fftPSF,x_im,y_im] = psfCamera1.ComputeFFTpsf(nSample,ntime);
% [fftPSF,x_im,y_im] = psfCamera1.ComputeFFTpsf();

%% PLOT THE PSF
plotType='surf';
limit=20*1e-3; % 
limit=[];

% plotType='contour';
% wave1=wave0;

%Select a wavelength
wave0=500; %nm

figure

subplot(1,2,1)
psfCamera1.drawFFTpsf(wave0,plotType,limit)
title (['Wavelength: ',num2str(wave0),' [nm]', ' NOT IN FOCUS'])

% SET the film in focus for the given wavelength

psfCamera1.autofocus(wave0,'nm',n_ob,n_im);
% Then Updata the 'Black Box Model'
bbmNEW=psfCamera1.bbmCreate(n_ob,n_im);

defocusNEW=psfCamera1.get('bbm','defocus','#wave');
paNEW=psfCamera1.get('bbm','primaryaberration','#wave');
%Re-Compute PSF
% ntime=8;  % window width= ntime*normalized_ExPdiameter
[fftPSF,x_im,y_im] = psfCamera1.ComputeFFTpsf(nSample,ntime);

%%PLOT
subplot(1,2,2)
psfCamera1.drawFFTpsf(wave0,plotType,limit)
% psfCamera1.drawFFTpsf(wave1,plotType)
title (['Wavelength: ',num2str(wave0),' [nm]', ' IN FOCUS'])

% % % set text box
textF='true';
imagePoint  = psfCamera1.get('bbm','gaussian image point')
figure
lens.draw
% psfCamera1.draw

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


