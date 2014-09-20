%Compute PSF for the given Pupil Function- Monochromatic transformation

%NOTE: With FFT2 method is the aliasing. In order to prevent this, the
%pupil must be small in comparison with the sampled area  In order to prevent this, the
% pupil must be small in comparison with the sampled area) 
% This is the common approach used in commercial ray tracing packages, such as Zemax or
% Code V


function [PSF]=psfPupil2PSF(Pupil,varargin)

%%INPUT
%Pupil: vector or matrix  for pupil function
%varargin {1}: type of illumination ['coherent','incoherent']
%               by default :'incoherent'

%OUTPUT
%PSF: psf [vector for vector pupil] [matrix for matrix pupil]

%% CHECK INPUT

if nargin>1
    ill_type=varargin{1};
else
    ill_type='incoherent';
end

%Scalar o vector psf?

if (ndims(Pupil)==2) 
    if size(Pupil,1)>1 && size(Pupil,2)>1
        fft_type='2d';
    else
        fft_type='1d';
    end
else
    error(['Not accepted matrix of ',num2str(ndims(Pupil)),' dimension'])
end

%% COMPUTE PSF

switch fft_type
    case {'1d'}
        PSF0=fftshift(fft(Pupil));
%         PSF0=(fft(Pupil));
    case {'2d'}
        PSF0=fftshift(fft2(Pupil));
%         PSF0=(fft2(Pupil));
end

%Illumination condition
switch ill_type
    case {'coherent';'coher'}
        PSF = PSF0;
    case{'incoherent';'incoher'}
        inten = (PSF0 .* conj(PSF0));   %intensity
        psf= real((inten));
        %DEBUG
        psf1=abs(PSF0).^2;
        %Normalize to have all sum to 1
        PSF = psf/sum(sum(psf));
    otherwise
        error([ill_type,' is not valid as illumination condition'])
end

% PSF0=(fft2(Pupil(:,:,li)));
% weight0=1./((wave(li,1).*defocusZ(li,1)).*sqrt(areaExP(li,1)));
% Coherent PSF
% PSF_coh(:,:,li)=weight0.*PSF0;
% Incoherent PSF
% inten = (PSF0 .* conj(PSF0));   %intensity
% psf= real(fftshift(inten));
% psf = psf/sum(sum(psf));
% PSF_in(:,:,li)=psf;
%     
% %     PSF_in(:,:,li)=weight0.^2*abs(PSF0).^2;
% end