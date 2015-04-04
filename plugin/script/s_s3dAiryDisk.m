%s_s3dAiryDisk.m
%
% Show a point source that is out of focus first.  Then move the film plane
% to the in focus position and compute the point focus again.  This will
% produce an Airy Disk.
%
% A lens is a set of surfaces without a point or film.  The properties of
% the lens independent of the point or film are computed with the lens
% function.  These include properties such as the ABCD matrix.  In the lens
% BBM, there is a relatively few number of properties that depend only on
% the lens itself. 
%  
% A psfCamera has the lens and a point and a film plane.  This permits MP
% to calculate many more aspects of the imaging system that require
% knowledge of the film plane and the object point.  For example, the
% Seidel aberrations for that point can be calculated.  So the BBM for the
% psfCamera is has more information in it than the BBM for the lens
% alone.
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

%% Initialize ISET

ieInit

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

%Specify which surface works as diaphragm
apertureIndex=6;
lens.set('aperture index', apertureIndex);

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
film = filmC('position', filmPosition, 'size', filmSize, 'wave', wave, 'resolution', resolution);


%% Create psfCamera structure (point source + lens + film [sensor])

% The film is not yet in the perfect focus plane.  It is just some film
% somewhere.
camera = psfCameraC('lens', lens, 'film', film, 'pointSource', pointSource);
camera.film.position

camera.bbmCreate;

%% Compute and plot the PSF

% Computing and plotting the Seidel wavefront aberrations can be based on
% the code inside ComputeFFTpsf().

% These sampling and space issues need to be simplified and clarified in
% general.  The current situation is a hack
nSample = 512; % power of 2 for Fast Fourier Transform
ntime   = 32;  % window width= ntime*normalized_ExPdiameter
[psf, x_im, y_im] = camera.ComputeFFTpsf(nSample,ntime); %#ok<NASGU,ASGLU>
% tmp = psf(:,:,1);
% vcNewGraphWin;
% surf(x_im(1,:),y_im(1,:),tmp)

% Select a wavelength
inFocusWave = 500; %nm
plotType='surf';
limit=[];
camera.drawFFTpsf(inFocusWave,plotType,limit);  % Should recompute the PSF
title (['Wavelength: ',num2str(inFocusWave),' [nm]', ' NOT IN FOCUS'])

%% Autofocus

% Place the film distance to be in focus for this wavelength
camera.autofocus(inFocusWave,'nm',n_ob,n_im);
imagePoint = camera.lens.findImagePoint(pointSource,n_ob,n_im);

% These should be the same
disp(imagePoint(1,3))
disp(camera.film.position(3))

% Rebuild the bbm
% AND Re-Compute PSF
% This shouldn't really be necessary.  When we reposition the film, all the
% old stuff should be cleared!
camera.bbmCreate;
[psf,x_im,y_im] = camera.ComputeFFTpsf(nSample,ntime);

% PLOT
camera.drawFFTpsf(inFocusWave,plotType,limit)
title (['Wavelength: ',num2str(inFocusWave),' [nm]', ' IN FOCUS'])

%% END
