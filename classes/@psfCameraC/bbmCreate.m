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
lens=obj.lens;
film=obj.film;
pSource=obj.pointSource;

if nargin>1
    n_ob=varargin{1}; %refractive index in object space
    n_im=varargin{2}; %n_im refractive index in image space
else
    n_ob=1; n_im=1;
end

%% Equivalent Black Box Model of lens
%     OptSys = paraxAnalyze(lens0);
%  lens0.set('black box model',OptSys);    % Equivalent Black Box Model of lens

[OptSys]=lens.bbmCreate(n_ob,n_im); 
unit=paraxGet(OptSys,'unit');
% Build Imaging System composed by Optical System+Film+ pSource
lV=paraxGet(OptSys,'lastvertex'); % last vertex of the optical system
F.z=film.position(3)+lV;
% F.z=film0.position(3)+OptSys.cardPoints.lastVertex;
F.res=film.resolution(1:2);F.pp=film.size; %um x um

[ImagSyst]=paraxOpt2Imag(OptSys,F,pSource,unit); 

[ps_heigth,ps_angle,ps_zpos]=coordCart2Polar3D(pSource(1),pSource(2),pSource(3)); %get image coordinate in polar coordinate

pSpolar(1)=ps_heigth;pSpolar(2)=ps_angle;pSpolar(3)=ps_zpos;
%  Equivalent Black Box Model of Imaging System
obj.set('black box model',ImagSyst,pSpolar);

%%    SET OUTPUT/s
if nargout>0
    varargout{1}=ImagSyst;
end



