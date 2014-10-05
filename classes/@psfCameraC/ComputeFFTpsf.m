function [PSF,x_im,y_im] = ComputeFFTpsf(obj,varargin)
% Compute the PSF using the Fourier Optics method: PSF=|PupilFunction|^2
%
% If you specify the number of samples, select a power of 2 (for Fast
% Fourier Trasform computation)  
%
% ASSUMED no APODIZATION (pupil illumination is uniform)
%
%   psfCamera.ComputeFFTpsf(nSample, ntime,unit)
%
%Examples
%   1:    ComputeFFTpsf(256,4,'mm')
%   2:    ComputeFFTpsf(256,4)
%   3:    ComputeFFTpsf(256)
%   4:    ComputeFFTpsf()
%
% The camera has a point source, lens, and film.
%
% INPUT
% varargin  {1}: number of sample of pupil function  [128x128 by default] [
% varargin  {2}:  % width of the window to sampling the pupil function  (normalized to ExP radius)  [from 1 to N] [2times by default]
% varargin  {3}: distance unit  ['mm' by default]
%
% OUTPUT
% varargout {1}:PSF = (nSample x nSample x num_wavelength)
% varargout {2}: x-coord (image space)  [Vector beacuse used equal sampling in both coordinate axis]
% varargout {3}: y-coord (image space) [Vector beacuse used equal  sampling in both coordinate axis]
%
% MP Vistasoft Team, Copyright 2014

%% TODO
% include a varargin{} to select type of apodization (or compute 
% illumination in real scene, i.e. partial occlusion of the PSF.
%
%

%% GET wavelength vector
wave=obj.get('wave'); % in nm
nW=size(wave(:),1); %number of sample
unit='mm';
wave_mm=wave*1e-6; % in mm
%% GET INPUT

% These are the default number of samples in one dimension of the pupil
% aperture. 
% def_ntim is the window is 2x the size of pupil.
% More comments to follow (MP)
def_nSample=128; def_ntim=2; def_unit='mm';
switch nargin
    case 1
        nSample=def_nSample;ntime=def_ntim;unit=def_unit;
    case 2
        nSample=varargin{1};
        ntime=def_ntim;unit=def_unit;
    case 3
        nSample=varargin{1};ntime=varargin{2};
        unit=def_unit;
     case 4
        nSample=varargin{1};ntime=varargin{2}; unit=varargin{3};  
    otherwise
        error('Too mucht inputs!!')        
        nSample=varargin{1};ntime=varargin{2}; unit=varargin{3};  
end


%% GET Coeffs for Pupils Function Polinomial
CoeffDef=obj.get('bbm','defocus'); %for defocus
CoeffPA=obj.get('bbm','primaryaberration'); %for primary aberration

% dZ_def=obj.get('bbm','defocusshift');
% NA=obj.get('bbm','numerical aperture');
% spot_rad_norm=dZ_def.*tan(NA); % radius for defocus


% Coeff1=0.5.*(NA./n_im).^2.*dZ; %Unique coeff

%% BUILD PUPIL FUNCTIONs
debugF='no-debug'; % flag for debugging

% I) Sampling space
range_pupil=2*ntime ; % al range od pupil function in normalized coordinate
d_pupil=range_pupil/nSample; %sampli
xn=[-nSample/2:nSample/2-1]*d_pupil;
yn=[-nSample/2:nSample/2-1]*d_pupil;
[Xn,Yn]=meshgrid(xn,yn);

Mask=circMask(sqrt(Xn.*Xn+Yn.*Yn)); 

%PUPIL NORMALIZED COORDINATE
% Vector
ro_v=sqrt(xn.*xn+yn.*yn);
theta_v= atan2(yn,xn);
%Matrix 
ro=sqrt(Xn.*Xn+Yn.*Yn); % radial coordinate [normalized to the Exit Pupil radius]
theta=atan2(Yn,Xn); % theta  

% II) Apodization (assumed uniform for all the wavelength)
typeApod='uniform';
[ApodW]=psfGetApodFun(typeApod,ro,theta,0);
% vcNewGraphWin; imagesc(ApodW)

% III) Phase Function = Seidel Aberration + Defocus
% Wavelength dependence (included in the function)
W_def=psfGetPhaseFun(CoeffDef,'defocus',ro,theta);
W_pa=psfGetPhaseFun(CoeffPA,'primary aberration',ro,theta);

% FOR DEBUG:
% W_def=psfGetPhaseFun(CoeffDef,'defocus',ro,theta,'debug',1);
% W_pa=psfGetPhaseFun(CoeffPA,'primary aberration',ro,theta,'debug',1);
% %how many wavelength of aberration
% CoeffDef_wave=paUnit2NumWave(CoeffDef);
% CoeffPA_wave=paUnit2NumWave(CoeffPA);

% IV) PUPIL FUNCTION= APODIZATION * EXP (-i*k*PHASE_FUNCTION)
% Wavelength dependence (included in the function)
[PupilFun]=psfGetPupilFun(ApodW,W_pa,W_def,wave_mm);
% vcNewGraphWin; imagesc(abs(PupilFun(:,:,1)))

%% COMPUTE POINT SPREAD FUNCTION (PSF)
light_type = 'incoherent';
[PSF] = psfPupil2PSF(PupilFun,light_type);

%% Get GRID sampling for Image Space
% Get ImagSyst
% ImagSyst=obj.get('imagingsystem');
% % Get point Source
% pSource=obj.get('pointsource');
% % Create a dummy object to compute the Numerical Apertuee
% Obj.z=pSource(3)+paraxGet(ImagSyst,'lastVertex'); 
% Obj.y=sqrt(pSource(1).^+pSource(2).^2); % eccentricity (height)
% % Numerical Aperture
% [NA]=paraxEstimateNumAperture(ImagSyst,Obj);
NA = obj.get('bbm','numerical aperture');
% Sampling GRID for image Space
coordType = 'normalized polar';

[x_im,y_im] = psfEstimateImageCoord(coordType,ro_v,theta_v,wave_mm,NA);


%% SET OUTPUT

obj.set('fft psf modulus',PSF);
coord.x = x_im;
coord.y = y_im;
obj.set('fft psf coordinate',coord);

% if nargout>1
%     varargout{1} = PSF;
%     varargout{2} = x_im;
%     varargout{3} = y_im;
% end


end