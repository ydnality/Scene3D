%% COMPARE THE POINT SPREAD FUNCTION computed by
% RayTracing Method + simulated diffraction(Andy's method) and 
% Wavefront estimation through paraxial approx (Michael's method)

% MISSING THE TELEPHOTO LENS BELOW.  Not sure where it went.  Ask MP.

%% INITIALIZATION
clear all
close all
ieInit

%% EVALUATION CONDITION
 
%Because I compute the PSF with the approx of the diffraction integral as
%the FFT, there are same constrain to be satisfy:
%1): number of sample has to be a power of 2
nSample=2^8; % in each direction, so the PSF for this value will 128x128 sample

%2)The pupil function has to be sample with a bigger window (at least 2
%times bigger). The following parameter specify how many time this window is considered
ntime=8;   

%IDEAL CASE  nSample->Inf   ntime->Inf, a trade-off is necessary

%% EXPERIMENTAL CONDITIONs
unit='mm';
wave=[400:50:600]*1e-6; % in mm 

%Select object and image space medium (air by default)
[n_ob]=RefractiveIndexDispersion(wave,unit,'air0');
[n_im]=RefractiveIndexDispersion(wave,unit,'air0');

%Select Optical System (2 Choices)
file_name='F_4_Telephoto Camera Lens'; % more info at pag 549 of the  attached BOOK
% file_name='F3_5_Cooke_Triplet.mat'; % % more info at pag 118 and 119 of the  attached BOOK

%Aperture (leave empty to use the default ones) 
% SUGGESTION if you want to change that, be careful to set low value otherwise you
% could enlarge the aperture and a lens will play the role of the
% diaphraghm
Aper_Diam=[]; %in mm

%Select the film
%position (any optical system has a default position of the film, but you
%can place the film in focus for a given wavelenght (@ 550 nm )
filmParam.waveRef=550*1e-6;
filmParam.mode='default';
% filmParam.mode='in-focus';
filmParam.dim=[36,48]; %Dimension: 36mm and 47mm along the two direction
filmParam.pitch=[0.01,0.01]; % pixel distance along the two directionum x um
%% SOURCE CONDITIONs

Obj.z=-10000; %object distance (has to be negative) [mm]
%you can specify the eccentricity as 
Obj.y=0; %height
%just leave one of the two field empty


%% STEP1: CREATE AN  IMAGING SYSTEM according to MICHAEL's standard
[ImagSyst]=paraxGenerateImagSystem(file_name,Aper_Diam,Obj,filmParam,n_ob,n_im,wave,unit);





%% STEP4:  ESTIMATE POINT SPREAD FUNCTION from ImagSyst (Michael's Method)

[PSF,x_im,y_im,Pupil,ExP_Diam,NA]=psfEstimate(ImagSyst,nSample,ntime); %Notice: the information of the object and film is included in the ImagSy


    %% SEVERAL PLOT
%select wavelenght to plot from the variable 'wave'
wave0=550*1e-6; %550 nm
inW0=find(wave==wave0);

%Plot of the PUPIL FUNCTION

%this is less relevant for you, but describe the wavefront error at the
%Exit Pupil (the deviation of the Optical Path of the Ray to the reference
%sphere, which is centered in the Gaussian image point and its radius is
%the distance between the Gaussian image point and the ExP center

%Normalized Pupil coordinate
range_pupil=2*ntime ; % al range od pupil function in normalized coordinate
d_pupil=range_pupil/nSample; 
xn=[-nSample/2:nSample/2-1]*d_pupil;
yn=[-nSample/2:nSample/2-1]*d_pupil;

lim1=ExP_Diam(inW0);
limit1=[lim1,lim1];
h1=figure(1)
[out1]=psfPlotPUPIL(Pupil,xn,yn,ExP_Diam,limit1,'surf',wave,wave0,h1);

%Plot PSF
limit2=[0.1,0.1]; %limit in [mm] of plot [scalar:limit equal in both direction, 1x2 to specify each direction

h2=figure(2)
out2=psfPLOT(PSF,x_im,y_im,limit2,'surf',wave,wave0,h2);

% Get a 1D profile for the specific wavelength
inW0=find(wave==wave0);
[Psf1D,vettX]=psfGet1Profile(PSF,x_im,y_im,wave,wave0,'peak-x');

figure(3)
plot(vettX,Psf1D)
title('PSF profile along X-axis')
xlabel('x [mm]'),ylabel('Normalized intensity')


% [Xim,Yim]=meshgrid(x_im(1,:),y_im(1,:))



%% STEP 5: ESTIMATE THE GEOMETRICAL PSF and with SIMULATED diffraction
% TO be completed by ANDY

%suggestion A: see within a loop for how many ray are need to come closer to
%the solution I found

%suggestion B: you can exploit the knowledge of the following parameter
% ANY POSITION along the OPTICAL AXIS is CONVERTED TO SUIT WITH YOUR COORDINATE SYSTEM
% (zero position: for the closer refractive surface to the film,positive sign in the way of the film)

z0=ImagSyst.surfs.list{ImagSyst.surfs.order(end)}.z_pos; % Your ZERO POSITION

%PARAMETER 1: Entrance Pupil
EnPob=ImagSyst.object{end}.Radiance.EnP; %the object
EnP_z=EnPob.z_pos-z0; %EnP position (it's a column vector becuase its wavelength dependent)
EnP_diam=EnPob.diam(:,1)-EnPob.diam(:,2); % EnP diameter (it's a column vector becuase its wavelength dependent)
%you could select the one wavelenght (or which show the bigger diameter)
%and use this information to shoot less rays
%EXAMPLE
% EnP_2angle=EnP_diam./EnP_z;
% [aM,inM]=max(abs(EnP_2angle));
% EnP_z=EnP_z(inM,1);
% EnP_diam=EnP_diam(inM,1);

%PARAMETER 2: Position of the film for 'BEST' FOCUS
bestFocus_z=ImagSyst.object{end}.ConjGauss.z_im-z0;
%its wavelenght dependent (select one)


%PARAMETER 3: Numerical Apertur
% it is more accurate to find DIFFRACTION LIMIT
NumApert=NA; % 
diffLimit_Radius= 1.22*wave'./NA; %it'is wavelenght dependent

%NOTE: the numerical aperture is defined as  NA=n_im*sin (alfa), where alfa is
%the angle subtended by the ExP radius for a sphere center in the
%Gaussian image point. In other word  alfa=tan(ExP_Diam/2*1/(distance
%ExP-GaussianImage point) or alfa=atan(effectiveF-number/2)
% so the NA=1/(2*effFnum) only for paraxial approximation (tan(k)=k,
% sin(k)=k) and the image space is in AIR


%% STEP2: CREATE AN IDENTICAL IMAGE SYSTEM according to ANDY's standard
[A_surfaceArray,A_lens,A_film,A_object]=paraxCreateScene3DSystem(ImagSyst,Aper_Diam);

%manually set film to focal position because the calculated one isn't
%working for some reason

bestFocus_z=ImagSyst.object{end}.ConjGauss.z_im-z0;

A_film.position = [0 -1.3 bestFocus_z(4)];
A_film.size = [.6 .6];
%A_film.wave = [550];

A_lens.diffractionEnabled = true;
%A_lens.set('wave', [550]);
% %DEBUG
psfCamera = psfCameraC('lens', A_lens, 'film', A_film, 'pointsource', A_object);
psfCamera.estimatePSF(100,true);
oi = psfCamera.oiCreate;
vcAddObject(oi); oiWindow; 


% ANYLIZE psf
yRange=psfCamera.film.size(1); xRange=psfCamera.film.size(2);
ny=psfCamera.film.resolution(1);nx=psfCamera.film.resolution(2);
nwave=psfCamera.film.resolution(3);


yv=[-ny/2:ny/2-1]*yRange*1e-3; xv=[-nx/2:nx/2-1]*xRange*1e-3;
yV=repmat(yv,nwave,1);xV=repmat(xv,nwave,1);

PSFandy=psfCamera.film.image;
h9=figure(9)
out9=psfPLOT(PSFandy,xV,yV,limit2,'surf',wave,wave0,h9);

% Get a 1D profile for the specific wavelength
inW0=find(wave==wave0);
[Psf1Dandy,vettXandy]=psfGet1Profile(PSFandy,xV,yV,wave,wave0,'peak-x');

figure(10)
plot(vettXandy,Psf1Dandy)
title('PSF profile along X-axis')
xlabel('x [mm]'),ylabel('Normalized intensity')


figure(11)
plot(vettX,Psf1D/max(Psf1D)),hold on
plot(vettXandy,Psf1Dandy/max(Psf1Dandy),'r')
title('PSF profile along X-axis')
xlabel('x [mm]'),ylabel('Normalized intensity')
