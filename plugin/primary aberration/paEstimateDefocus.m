function [defocusCoeff,varargout] = paEstimateDefocus(ImagSyst,Obj,varargin)
%Compute the Defocus Coeffs for the normalized coordinate
% based on different level of approxiamtion of Nijboer-Zernike Theory
% ADVICE: select 'highNA' as method to compute NA
%
% function [defocusCoeff] = paEstimateDefocus(ImagSyst,Obj,varargin)
%
%INPUT
%  ImagSyst: image system struct
%  Obj:      point source object structure  {.z: position along optical axis; .y:height eccentricity}
%  vararging {1}: method of approximation of the defocus phase aberration
%                  function
%
%OUTPUT
%   Coeff: coeffs for defocus (accordin to the selected method can be from 1
%          to 3 coeffs)
%    varargout{1}: distance between gaussian image point and film position
%
% MP Vistasoft 2014


%% See APPENDIX A: 
%%Braat, Joseph, Peter Dirksen, and Augustus JEM Janssen. "Assessment of an extended Nijboer–Zernike approach for the computation of optical point-spread functions." JOSA A 19.5 (2002): 858-870.
% Wf=dZ(1-sqrt(1-sin(p)^2); p=ro*NA
% 
%and 'The Extended Nijboer-Zernike Driffraction Theory and its
%Applications' Sven van Haver Thesis

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


%% Get parameters from the Imaging System
if size(ImagSyst.film{end}.z_pos,1)==nW
%     defocusZ(:,1)=ImagSyst.film{end}.z_pos; %last film
    defocusZ(:,1)=paraxGet(ImagSyst,'film position');  %last film
else
    defocusZ(:,1)=repmat(paraxGet(ImagSyst,'film position'),[nW,1]);
end
if isinf(ImagSyst.object{end}.z_pos)
    defocusZ(:,2)=paraxGet(ImagSyst,'image focal point'); 
%     defocusZ(:,2)=ImagSyst.cardPoints.dFi+ImagSyst.cardPoints.lastVertex; 
else
    defocusZ(:,2)=paraxGet(ImagSyst,'imagepointposition');
%     defocusZ(:,2)=ImagSyst.object{end}.ConjGauss.z_im;
end


for li=1:nW
    ExPDiam(li,1)=ImagSyst.object{end}.Radiance.ExP.diam(li,1)-ImagSyst.object{end}.Radiance.ExP.diam(li,2);
    fi=paraxGet(ImagSyst,'effectivefocallength');
%     efl(li,1)=ImagSyst.cardPoints.fi(li,1);
    efl(li,1)=fi(li,1);
    [NA(li,1)]=paraxNumAperture(ExPDiam(li,1),efl(li,1),ImagSyst.n_im(li,1));
     if nargin >3
         approx_method=varargin{1};
     elseif NA(li,1)<=0.6
         approx_method='1term';
     else
        approx_method='2term';
        % approx_method='3term';
     end
    [Coeff(li,:),dZ(li)]=paEstimateDefocusCoeff(defocusZ(li,1),defocusZ(li,2),NA(li,:),ImagSyst.n_im(li,:),approx_method);
end


%% SET OUTPUT
defocusCoeff.wave=ImagSyst.wave;
defocusCoeff.unit=ImagSyst.unit;
nterms=size(Coeff,2); % number of terms
defocusCoeff.type='defocus';
defocusCoeff.note=[num2str(nterms),'term'];
for ni=1:nterms
    field_name=['W',num2str(ni*2),'0'];
    defocusCoeff=setfield(defocusCoeff,field_name,Coeff(:,ni)); % set defocus term
end

if nargout>1
    varargout{1}=dZ;
end
