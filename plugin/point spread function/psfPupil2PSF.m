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
    waveDep=0; %wavelength-dependence [FALSE]
else
    fft_type='2d';
    waveDep='true'; %wavelength-dependence [TRUE]
%     error(['Not accepted matrix of ',num2str(ndims(Pupil)),' dimension'])
end

%% COMPUTE PSF


if waveDep
    nW=size(Pupil,3); % number of wavelength
    
    for li=1:nW       
     switch fft_type
            case {'1d'}
                PSF0(:,:,li)=fftshift(fft(Pupil(:,:,li)));
        %         PSF0=(fft(Pupil));
            case {'2d'}
                PSF0(:,:,li)=fftshift(fft2(Pupil(:,:,li)));
        %         PSF0=(fft2(Pupil));
        end

        %Illumination condition
        switch ill_type
            case {'coherent';'coher'}
                PSF(:,:,li) = PSF0(:,:,li);
            case{'incoherent';'incoher'}
                inten = (PSF0(:,:,li) .* conj(PSF0(:,:,li)));   %intensity
                psf= real((inten));
                %DEBUG
                psf1=abs(PSF0(:,:,li)).^2;
                %Normalize to have all sum to 1
                PSF(:,:,li) = psf./sum(sum(psf));
            otherwise
                error([ill_type,' is not valid as illumination condition'])
        end
    end
    
else
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
end


