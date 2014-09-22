% Compute the PSF for the given zernike/4th order coeffs, specifying the
% ExP radius, number of sampling, the defocus (if present) and eventually
% apodization

% 
function [PSF_in]=psfCompute(AberCoeff,ExP_Diam,nSample,timesP,defocusZ,varargin)

%INPUT
%AberCoeff: struct with Zernike or PeakCoeff
%ExP_Diam: Exit pupil Diameter[scalar or wavelenght dependent{Nx1}]
%nSample: number of sample for axis( scalar: uniform samplein, 2x1 specify [numXsample,numYsample]
%timesP: n-times the ExP diameter for sampling the Pupil function
%defocusZ:[Nx2] [distance between -ExP and the Gaussian Image point, distance
%between -ExP and the SENSOR]            SENSOR ASSUMED FLAT
%Varargin  {1}: apodization type
%           {2}:  if apodization is 'customized' specify the parameter

%OUTPUT
%PSF: wavelength dependent point spread function for incohrent illumination

%NOTE:
% N: number of wavelength sample
%Assumed CIRCULAR APERTURE


%% Obsolete function: substitued by the function ''

%% GET FONDAMENTAL PARAMETERs
wave=AberCoeff.wave;
nw=size(wave,2); %number of sample
unit=AberCoeff.unit;

if nargin>4
    apod_type=varargin{1};
else
    apod_type='uniform';
end


%% CREATE SAMPLE along the axis [normalized] xn[-1 0 1]; yn[-1 0 1];
if length(nSample)==1
    numx=nSample; numy=nSample;
else
    numx=nSample(1);numy=nSample(2);
end
dx=1./(numx/2-1);
dy=1./(numy/2-1);

%normalized axis
Rmax=timesP; %normalized radius 1
xn_pos=[0:dx:Rmax]; xn_neg=[-Rmax:dx:-dx]; xn=[xn_neg,xn_pos];
yn_pos=[0:dx:Rmax]; yn_neg=[-Rmax:dx:-dx]; yn=[yn_neg,yn_pos];

[Xn,Yn]=meshgrid(xn,yn);

%% CREATE MASK
%Assumed uniform
Mask=circMask(sqrt(Xn.*Xn+Yn.*Yn));

%% CREATE APODIZATION FUNCTION
switch apod_type
    case {'uniform'}
        Apod=circMask(sqrt(Xn.*Xn+Yn.*Yn));
    case {'Stiles-Crawford'}
        error(['The method ',apod_type,' for apodization has to be completed!'])
    case {'Gaussian';'gaussian'}
        if nargin>5
            sigma=varargin{2};
        else
            sigma=1; %by default sigma is equal to ExP Radius
        end
        Apod=fspecial('gaussian',[length(xn), length(yn)],sigma);
    case {'customized'}
        error(['The method ',apod_type,' for apodization has to be completed!'])
    otherwise
end

%% CREATE PHASE WAVEFRONT ABERRATION
fname=paGetFieldNames(AberCoeff);
numF=length(fname); %number of coeffs
%check is available Zernike or Peak Coeffs
switch fname{1}(1)
    case {'C','c'}
        wvf_type='zernike';
    case {'W','w','A','a'}
        if length(fname{1})==3
            wvf_type='peak';
        else
            error([fname{1},' is not a valid type of aberration. Probabily the peak value has to be computed '])
        end
    otherwise
        error([fname{1},' is not a valid type of aberration '])
end

switch wvf_type
    case {'zernike'}
        
    case {'peak'}
        ro=sqrt(Xn.*Xn+Yn.*Yn); %normalized radius
        theta=atan2(Yn,Xn); % theta  
%         theta=atan(Yn./Xn); % theta  
        for fi=1:numF
            kC=str2num(fname{fi}(2)); %get radius order
            lC=str2num(fname{fi}(3)); %get cosine-theta order
            
            Wkl(:,fi)=getfield(AberCoeff,fname{fi});%weight 
%             B(:,:,fi)=pa4thWaveAber(kC,lC,ro,theta); %bases 
            B(:,:,fi)=Mask.*pa4thWaveAber(kC,lC,ro,theta); %bases 
%             %DEBUG
%             figure(fi)
%             surf(B(:,:,fi))
%             title(fname{fi})
            for li=1:nw
                WB(:,:,fi,li)=Wkl(li,fi).*B(:,:,fi);
            end
        end
        PhaseW=sum(WB,3);
        %
end
        
        

%% DEFOCUS TERM
% used defocused pupil function described by Mahajan- PartII Wave
% Diffraction Optics, equation 1-47

coeffDefocus=0.5*(1./defocusZ(:,1)-1./defocusZ(:,2)).*(ExP_Diam/2).^2;
baseDefocus=pa4thWaveAber(2,0,ro,theta); %bases
for li=1:nw
    PhaseDefocus(:,:,li)=coeffDefocus(li,1).*baseDefocus;
end


%% PUPIL FUNCTION

for li=1:nw
    Kl(li)=2*pi./wave(li,1);
%     Pupil(:,:,li)=Mask.*Apod.*exp(-i*Kl(li).*(PhaseW(:,:,li)+PhaseDefocus(:,:,li)));
    Pupil(:,:,li)=Apod.*exp(-i*Kl(li).*(PhaseW(:,:,li)+PhaseDefocus(:,:,li)));
end



%% PSF
areaExP=pi*(ExP_Diam./2).^2;
for li=1:nw
    PSF0=(fft2(Pupil(:,:,li)));
    weight0=1./((wave(li,1).*defocusZ(li,1)).*sqrt(areaExP(li,1)));
    %Coherent PSF
    PSF_coh(:,:,li)=weight0.*PSF0;
    %Incoherent PSF
    inten = (PSF0 .* conj(PSF0));   %intensity
    psf= real(fftshift(inten));
    psf = psf/sum(sum(psf));
    PSF_in(:,:,li)=psf;
    
%     PSF_in(:,:,li)=weight0.^2*abs(PSF0).^2;
end
