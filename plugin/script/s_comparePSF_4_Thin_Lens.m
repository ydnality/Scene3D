% s_comparePSF_4_Thin_Lens.m
%
%  Check the diffraction limited calculation for a thin lens and a point
%  source based on MP's black box methods
%
% MP, Vistasoft Team, 2015

%% INITIALIZATION
clear all
close all
ieInit

%% SET
wave = (400:50:600)'*1e-6; %mm
nWave   = size(wave,1);
unit = 'mm';

%Refractive index object space
n_ob=ones(nWave,1);
n_im=ones(nWave,1);

%THIN LENS PARAMETER 
%
thinLens.z_pos = 0; % the thin lens is in the reference position
thinLens.diam  = 3; %diameter of thin lens (which play as diaphragh 

efl0 = 17.1;       % mm effective focal length (like human eye)
thinLens.efl = repmat(efl0,size(wave,1),1); %ASSUMED NO DISPERSION

% F-number (adequate for diffraction in paraxial approximation)
thinLens.effFnum=thinLens.efl./thinLens.diam; %effective F-number

%NA: numerical aperture (more useful in general case of diffraction)
[thinLens.NA] = paraxNumAperture(thinLens.diam,thinLens.efl,n_im);

% Diffraction limit (More accurate if computed with NUMERICAL APERTURE
thinLens.diffLimit_Radius=1.22*(wave./(2*thinLens.NA));
% thinLens.diffLimit_Radius=1.22*(wave.*thinLens.effFnum);%

%% Film position
%case:in focus
film_z=efl0;
%case:not in focus
def_Z=0; %mm
film_z=efl0+def_Z;


%% MICHAEL's APPROACH

    %%STEP 1: create a pupil function
    
    %Because I compute the PSF with the approx of the diffraction integral as
    %the FFT, there are same constrain to be satisfy:
    %1): number of sample has to be a power of 2
    nSample=2^8; % in each direction, so the PSF for this value will 128x128 sample

    %2)The pupil function has to be sample with a bigger window (at least 2
    %times bigger). The following parameter specify how many time this window is considered
    ntime=8;   

    %IDEAL CASE  nSample->Inf   ntime->Inf, a trade-off is necessary

%Normalized coordinate
range_pupil=2*ntime ; % in normalized coordinate
d_pupil=range_pupil/nSample;
xn=[-nSample/2:nSample/2-1]*d_pupil;
yn=[-nSample/2:nSample/2-1]*d_pupil;


%Related Image Space normalized coordinate
d_image=1/range_pupil;
image_range=1/d_pupil;
xim=[-nSample/2:nSample/2-1]*d_image;
yim=[-nSample/2:nSample/2-1]*d_image;



    %% STEP 2 : Generate a  WAVEFRONT for the thin lens
list_wave={'W040';'W131';'W222';'W220';'W311'};
% The aberrazion are expressed as number of wave
num_wave0=0; %aberration free 
num_wave=num_wave0*[1,1,1,1,1]; %NO ABERRATION for the thin lens


    %% CREATE  a circular mask to the normalized pupil coordinate 
Mask=circMaskXY(xn,yn);
    
    %% CREATE  the Apodization Factor (intensity of the ray at the pupil plane) m
    % Assumed UNIFORM
Apod=circMaskXY(xn,yn);

    %% SET ABERRATION COEFFS and CREATE WAVEFRONT ERROR FUNCTION
    WaveCoeff_wave=struct;
for ci=1:length(list_wave)
    WaveCoeff_wave=setfield(WaveCoeff_wave,list_wave{ci},num_wave(ci));
end
WaveCoeff_wave.wave=wave;
WaveCoeff_wave.unit=unit;
%Convert wavefront from #wave to unit
[WaveCoeff]=paNumWave2Unit(WaveCoeff_wave);

%Compute Peak Value Coeff
ymax=1;y=1;
[PeakCoeff]=paWave4thOrder2PeakValue(WaveCoeff,y,ymax);

[PhaseW,ro]=paComputePhaseWvf(PeakCoeff,xn,yn,Mask);     
    
    %% STEP 3: Add a DEFOCUS TERM


if size(film_z,1)==nWave
    defocusZ(:,1)=film_z;defocusZ(:,2)=thinLens.efl;    
else 
    defocusZ(:,1)=repmat(film_z,nWave,1);
    defocusZ(:,2)=thinLens.efl;    
end
%Compute Wavefron term for defocus
for li=1:nWave
    PhaseDefocus(:,:,li)=paDefocus4thWaveAber(defocusZ(li,1),defocusZ(li,2),thinLens.NA (li,:),ro,'high');
end
   
    
   
    %% STEP4: ASSEMBLE ALL THE TERMs
    
    



for li=1:nWave
    %Phase scaling
    Kwave=-i*2*pi/wave(li);
    PhW0=Apod.*exp(Kwave.*(PhaseW(:,:,li)+PhaseDefocus(:,:,li)));

    %% COMPUTE THE PSF
    PSF(:,:,li)=psfPupil2PSF(PhW0,'incoherent');
    %Find 'real' coord in the image plane
    [xIM(li,:),yIM(li,:)]=psfNormalized2RealCoordinate(xim,yim,wave(li),thinLens.NA(li));
end


%% PLOT THE RESULTING PSF

%Select witch wavelength plot
wave0=550*1e-6; %550 nm
limit1 = 2*[thinLens.diam,thinLens.diam];
h1 = vcNewGraphWin;
out1=psfPLOT(PSF,xIM,xIM,limit1,'surf',wave,wave0,h1);


%%  PLOT 1D PROFILE of the PSF along the its Max VALUE

[Psf1D,vettX]=psfGet1Profile(PSF,xIM,yIM,wave,wave0,'peak-x');

% DIFFRACTION LIMIT by Theory
% v=r*lambda/NA
%I(v)/I0=[2*J1(v)/v]^2
%PSF_incoherent light=abs(I(v)/I0)^2;
inW0=find(wave==wave0);
v=xIM(inW0,:)*(thinLens.NA(inW0))/wave0;
Psf1D_DiffLim=abs(somb(2*v)).^2;
%  Psf1D_DiffLim=airy(xIM(inW0,:));

vcNewGraphWin; 
plot(vettX,Psf1D./max(Psf1D),'LineWidth',2)
hold all
plot(vettX,Psf1D_DiffLim./max(Psf1D_DiffLim),'--r','LineWidth',2)
title('PSF profile along X-axis')
xlabel('x [mm]'),ylabel('Normalized intensity')
legend ('PSF estimated by FFT of Pupil Function','Diffraction Limited PSF (Theoretical)')


%% PLEASE ANDY INSERT HERE  YOU ESTIMATION