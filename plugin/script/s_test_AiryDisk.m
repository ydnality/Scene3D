%s_test_Airy disk.m
%
%  A lens is a set of surfaces without a point or film.  The properties of
%  the lens independent of the point or film are computed with the lens
%  function.  These include properties such as the ABCD matrix.  Maybe some
%  other stuff.  In the lens BBM, there is a relatively few number of
%  properties that depend only on the lens itself.
%  
%  When you make a psfCamera, it has the lens above and also a point and a
%  film plane.  This permits MP to calculate many more aspects of the
%  imaging system that require knowledge of the film plane and the object
%  point.  For example, the Seidel aberrations for that point can be
%  calculated.  So the BBM for the psfCamera is bigger with more
%  information in it than the BBM for the lens alone.
%
% PROBLEMS
%  There are serious problems when the point source distance and film
%  distance parameters are out of any sensible range.  We should do some
%  checking or set up the script so that really bad stuff can't happen in
%  these cases.
%
%Compute the wave PSF for a psfCamera (in Scene3D)
%
% The implemented method is the following:
% 1) Estimate the wavefront aberration (Seidel aberration)
% 2) Exploit the Fourier Optics, for incoherent light, to compute the PSF 
%
% MP Vistasoft, 2014

%% Initialize Iset
s_initISET

%% Specify the point source

% The center of the first camera aperture is usually at [0 0 0].
% The object positions data are in -z.  We are using a right-handed
% coordinate system. 

% Angulare eccentricity (in degree)
angle_field = 2;        % Field height, basically
ps_dist     = 0.25*1e6; % Good distance

% Bad distance.  Aliasing.  Problem connected to the number of 
% ps_dist = 0.25*1e2; 

% Point source in object space  and refractive indices of both object and
% image space 
wave = 400:50:600;

% Create point source coordinate
ps_x        =  ps_dist.*tan(angle_field/180*pi);
pointSource = [ps_x 0 -ps_dist];  % A very distant point.

%% Create lens

%Refractive indices on either side of the lens

% Default is air (vacuum).  For the eyeball n_im might be 1.3
n_ob = 1 ; % object space
n_im = 1;  % image space

% lensFile = fullfile(s3dRootPath, 'data', 'lens', '2ElLens.mat');
lensFile = fullfile(s3dRootPath, 'data', 'lens', 'dgauss.50mm.mat');
import   = load(lensFile,'lens');
lens     = import.lens;
lens.set('wave', wave);
n = lens.get('nArray');

% If the aperture is very small then the lens will be approximately
% diffraction limited.  A small number is 0.5.  A bigger number is 5.
% We see the Airy disk when 0.5 and not when 5.
lens.apertureMiddleD = 3;   %mm, Big
vcNewGraphWin;
lens.draw

%% Declare film
% filmPosition = [0 0 51.2821	];  % Good for 2Elens
% filmPosition = [0 0 107]; 
filmPosition = [0 0 35.5];          % Good for dgauss.50mm.  True focal about 37.3mm
filmDiag   = 10;  % Millimeters
filmSize   = [filmDiag/sqrt(2) filmDiag/sqrt(2)];
resolution =  [300 300 length(wave)];
film = pbrtFilmC('position', filmPosition, 'size', filmSize, 'wave', wave, 'resolution', resolution);

%% Create psfCamera structure (point source + lens + film [sensor])

% At this point the film is not yet in the perfect focus plane.  It is just
% some film somewhere.
psfCamera1 = psfCameraC('lens', lens, 'film', film, 'pointSource', pointSource);
psfCamera1.film.position

psfCamera1.bbmCreate;

% THIS FUNCTION IS NOT WORKING CORRECTLY.  SO, BW IS SCREWING AROUND TRYING
% TO UNDERSTAND WHY (October, 2014).
% Problem seems to be when the distance is too small (see up above).
% 
% psfCamera1.autofocus(wave0,'nm',n_ob,n_im);
% psfCamera1.film.position(3) = -1*psfCamera1.film.position(3);
% 
% imagePoint = psfCamera1.lens.findImagePoint(pointSource,n_ob,n_im)

%% Not sure why these were here.
% bbmOLD     = psfCamera1.bbmCreate(n_ob,n_im);
% defocusOLD = psfCamera1.get('bbm','defocus','#wave');
% paOLD      = psfCamera1.get('bbm','primary aberration');

%% SET all the primary aberration to 0

% We should eliminate this.  MP used it to set the aberrations to zero, but
% these should always be derived, really.
% psfCamera1.bbmSetAberration('all zero');

%% Compute and plot the PSF

% Computing and plotting the Seidel wavefront aberrations can be based on
% the code inside ComputeFFTpsf().

% These sampling and space issues need to be simplified and clarified in
% general.  The current situation is a hack
nSample = 512; % power of 2 for Fast Fourier Transform
ntime   = 32;  % window width= ntime*normalized_ExPdiameter
[fftPSF, x_im, y_im] = psfCamera1.ComputeFFTpsf(nSample,ntime); %#ok<NASGU,ASGLU>
% tmp = fftPSF(:,:,1);
% vcNewGraphWin;
% surf(x_im(1,:),y_im(1,:),tmp)

plotType='surf';
limit=[];
% Select a wavelength
vcNewGraphWin;
wave0=500; %nm
psfCamera1.drawFFTpsf(wave0,plotType,limit);  % Should recompute the PSF
title (['Wavelength: ',num2str(wave0),' [nm]', ' NOT IN FOCUS'])

%% Autofocus

% Set the film position into the focus for the given wavelength and update
% the black box model.
% This should produce an Airy Disk when the aperture diameter is small.
psfCamera1.autofocus(wave0,'nm',n_ob,n_im);

bbmNEW     = psfCamera1.bbmCreate(n_ob,n_im);
defocusNEW = psfCamera1.get('bbm','defocus','#wave');
paNEW      = psfCamera1.get('bbm','primary aberration');

% % Set all aberration to zero
% psfCamera1.bbmSetAberration('all zero');

%Re-Compute PSF
% ntime=8;  
% window width= ntime*normalized_ExPdiameter
%
[fftPSF,x_im,y_im] = psfCamera1.ComputeFFTpsf(nSample,ntime);

% PLOT
vcNewGraphWin;
psfCamera1.drawFFTpsf(wave0,plotType,limit)
% psfCamera1.drawFFTpsf(wave1,plotType)
title (['Wavelength: ',num2str(wave0),' [nm]', ' IN FOCUS'])

%% END
