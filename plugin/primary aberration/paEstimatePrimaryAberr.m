%Compute the Primary Aberration Coeffs for the normalized coordinate in
%power series expansion
% W40 ro^4+W31 ro^3 cos(T)+W22 ro^2 cos(T)^2+ W20 ro^2 + W11 ro cos(T)
% All coeff include high field (H)
% W40= sherical aberration = W040 
% W31= coma                = W131 H
% W22= astigmatism         = W222 H^2
% W20= field curvature     = W220 H^2 
% W11= distortion          = W311 H^3

% reference -Chapter 5- "Introduction to Aberrations in Optical Imaging System-J. Sasian"

function [Coeff] = paEstimatePrimAberr(ImagSyst,Obj,varargin)

%INPUT
%ImagSyst: image system struct
%Obj:      point source object structure  {.z: position along optical axis; .y:height eccentricity}
%vararging {1}: method of approximation of the defocus phase aberration
%function


%OUTPUT
%Coeff: set of coeffs for Primary Aberration (Power series expasion for
%wavefront aberration)


%% More detail in reference "Introduction to Aberrations in Optical Imaging System-J. Sasian"

%% CHECK INPUT
if isempty(Obj)
    error('The Point Source structure is empty')
end

%% GET Useful parameter
unit=ImagSyst.unit;
wave=ImagSyst.wave;
nW=size(wave); % number of wavelength

%% Add the Point Source to the Imaging  system
profile1='point';
[source1]=paraxCreateObject(Obj.z,Obj.y,profile1,unit);

[ImagSyst]=paraxAddObject2ImagSyst(ImagSyst,source1);

%% SET OUTPUT

Coeff=ImagSyst.object{end}.Wavefront.PeakCoeff;
Coeff.type='primaryaberration';

% %% Get parameters from the Imaging System
% if size(ImagSyst.film{end}.z_pos,1)==nW
%     defocusZ(:,1)=ImagSyst.film{end}.z_pos; %last film
% else
%     defocusZ(:,1)=repmat(ImagSyst.film{end}.z_pos,nW,1);
% end
% if isinf(ImagSyst.object{end}.z_pos)
%     defocusZ(:,2)=ImagSyst.cardPoints.dFi+ImagSyst.cardPoints.lastVertex; 
% else
%     defocusZ(:,2)=ImagSyst.object{end}.ConjGauss.z_im;
% end
% 
% 
% for li=1:nW
%     ExPDiam(li,1)=ImagSyst.object{end}.Radiance.ExP.diam(li,1)-ImagSyst.object{end}.Radiance.ExP.diam(li,2);
%     efl(li,1)=ImagSyst.cardPoints.fi(li,1);
%     [NA(li,1)]=paraxNumAperture(ExPDiam(li,1),efl(li,1),ImagSyst.n_im(li,1));
%      if nargin >3
%          approx_method=varargin{1};
%      elseif NA(li,1)<=0.6
%          approx_method='1term';
%      else
%         % approx_method='2term';
%         % approx_method='3term';
%      end
%     Coeff(:,:,li)=paEstimateDefocusCoeff(defocusZ(li,1),defocusZ(li,2),NA(li,:),ImagSyst.n_im(li,:),approx_method);
% end