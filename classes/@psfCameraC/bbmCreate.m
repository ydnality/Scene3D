function [varargout] = bbmCreate(obj,varargin)
% Create the black box model of the psfCamera structure of SCENE3D to the Imaging  System structure to
% get its analysis through paraxial approximation (first order optics)
%
%  ImagSyst = psfCamera.bbmCreate (varargin)
%
%
%INPUT
%   obj: lens object of SCENE3D
%   varargin {1}: n_ob refractive index in object space
%   varargin {2}: n_im refractive index in image space
%OUTPUT
%   varargout{1}=ImagSyst: Optical System structure.
%
% MP Vistasoft 2014


 

%% Get inputs
% lens=obj.lens; %NOT NEEDED in this function
% film=obj.film;  %NOT NEEDED in this function
pSource=obj.pointSource;

%%  CHECK number of INPUTs and set refractive indices of the medium
if nargin>1
    n_ob=varargin{1}; %refractive index in object space
    n_im=varargin{2}; %n_im refractive index in image space
else
    n_ob=1; n_im=1;
end

% %% Get Optical System from the lens in the psfCamera [ Michael's script]
% OptSyst=obj.get('optical system',n_ob,n_im);
% 
% unit=paraxGet(OptSys,'unit');
% 
% % Build Imaging System composed by Optical System+Film+ pSource
% 
% lV=paraxGet(OptSys,'lastvertex'); % last vertex of the optical system
% F.z=film.position(3)+lV;
% 
% F.res=film.resolution(1:2);F.pp=film.size; %um x um
% 
% % %% Equivalent Black Box Model of lens
% % lens.bbmCreate(n_ob,n_im); 
% % 
% % %Get Optical System from Michael's script
% % [OptSys]=lens.get('optical system');
% % 
% % unit=paraxGet(OptSys,'unit');
% 
% [ImagSyst]=paraxOpt2Imag(OptSys,F,pSource,unit); 

%% GET (by compute) the IMAGING SYSTEM
ImagSyst=obj.get('imaging system',n_ob,n_im);



% OLD FUNCTIONs [ps_heigth,ps_angle,ps_zpos]=coordCart2Polar3D(pSource(1),pSource(2),pSource(3)); %get image coordinate in polar coordinate
[ps_heigth,ps_angle,ps_zpos]=cart2pol(pSource(1),pSource(2),pSource(3)); %get image coordinate in polar coordinate


pSpolar(1)=ps_heigth;pSpolar(2)=ps_angle;pSpolar(3)=ps_zpos;
%  Equivalent Black Box Model of Imaging System
obj.set('black box model',ImagSyst,pSpolar);

%%    SET OUTPUT/s
if nargout>0
    varargout{1}=obj.get('black box model','all');
end



