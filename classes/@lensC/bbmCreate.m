function [varargout] = bbmCreate(obj,varargin)
% Create the black box model of the lens structure of SCENE3D to the Optical System structure to
% get its analysis through paraxial approximation (first order optics)
%
%  BBM = lens.bbmCreate ()
%    or
%  BBM = lens.bbmCreate (n_ob,n_im)
%    or
%   psfCamera = lens.bbmCreate ('all',pointSource,film)
%       or
%   psfCamera = lens.bbmCreate ('all',pointSource,film,n_ob,n_im)
%
%INPUT
%   obj: lens object of SCENE3D
%   varargin {1}: n_ob refractive index in object space
%   varargin {2}: n_im refractive index in image space
%   
%   or
%   varargin {1}: 'all'
%   varargin {2}: pointSource (3 value for x,y,z position)
%   varargin {3}: film (struct with following fields:  .;.)
%
%OUTPUT
%   varargout{1}= Black Box Model or psfCamera (if varargin{1}='all')
%
% MP Vistasoft 2014


%% CHECK INPUT and BUILD OPTICAL SYSTEM
psfFlag=0; %no psfCamera object to compute
if nargin>1    
    switch varargin{1}
        case {'all'}
            psfFlag=1; %YES psfCamera object to compute
            lens0=obj;            
            pSource=varargin{2}; %get point source
            film0=varargin{3};
            if nargin >4
                n_ob=varargin{4};n_im=varargin{5};
            else
                n_ob=1; n_im=1;
            end
            OptSyst=obj.get('optical system',n_ob,n_im);
        otherwise                   
            n_ob=varargin{1}; n_im=varargin{2};
            OptSyst=obj.get('optical system',n_ob,n_im);
%             OptSyst=obj.bbmComputeOptSyst(n_ob,n_im);
    end
else            
    n_ob=1; n_im=1;
    OptSyst=obj.get('optical system',n_ob,n_im);
%     OptSyst=obj.bbmComputeOptSyst(n_ob,n_im);
end

%% Append Optical System field to the Black Box Model of the lens
obj.set('black box model',OptSyst);  

if psfFlag
    psfCamera0 = psfCameraC('lens', lens0, 'film', film0, 'pointSource', pSource);
    psfCamera0.bbmCreate();
    varargout{1} =psfCamera0;
else
    varargout{1} =obj.get('bbm','all');
end



%% TODO
%
% Put paraxAnalyze into the @lensC directory as a function called
% bbmCreate.  And then the syntax will be
%
%   lens.bbmCreate
%  
%   Inside of the @lensC directory the function bbmCreate(obj) does what is
%   in paraxAnalyze.
%   
%  


% %% GET RELEVANT PARAMETERs for the COMPUTATION
% 
% unit='mm'; %unit
% wave=obj.wave*1e-6; % in mm
% nw=length(wave); %num wavelength
% 
% nelem=length(obj.surfaceArray);
% 
% %Initialize some vector
% N=ones(nw,nelem); 
% 
% %Useful parameter
% inD=1;
% 
% %% Get the parameter to build the Optical System
% for ni=1:nelem
%     %Get the structure
%     S=obj.surfaceArray(ni);
%     if all(S.n==0)
%         surftype{ni}='diaphragm';  
%         if (S.apertureD==obj.apertureMiddleD) %Check if the aperture size is changed
%             Diam(ni)=S.apertureD; %aperture diameter
%         else
%             Diam(ni)=obj.apertureMiddleD; %set aperture change
%         end
%         %save indices of the aperture
%         indDiaph(inD)=ni;
%         inD=inD+1;
% 
%         if ni>1
%             N(:,ni)=N(:,ni-1); %refractive indices
%        end
%     else
%         surftype{ni}='refr';
%         N(:,ni)=S.n';           %refractive indices                
%         Diam(ni)=S.apertureD; %aperture diameter
%     end
%     Radius(ni)=S.sRadius; %radius of curvature
%     posZ(ni)=S.get('zintercept');
%     %% OTHER FIELDs
%     %asphericity (conical parameter)
%     
% end
% % Set new origin as the First Surface
% PosZ=posZ-posZ(1);
% 
% %% Build several surface
% for k=1:length(Radius)
%     R(k)=Radius(k);
%     z_pos(k)=PosZ(k);
%     n(:,k)=N(:,k);
%     diam(k)=Diam(k);
%     switch surftype{k}
%         case {'refr'}          
%             surf{k}=paraxCreateSurface(z_pos(k),diam(k),unit,wave,surftype{k},R(k),n(:,k));
%         case {'diaphragm','diaph'}
%             surf{k}=paraxCreateSurface(z_pos(k),diam(k),unit,wave,surftype{k});
%         case {'film'}
%     end
% end
% 
% 
% %% CREATE OPTICAL SYSTEM and (possibly) IMAGING SYSTEM and psfCamera 

% [OptSyst] = paraxCreateOptSyst(surf,n_ob,n_im,unit,wave);


%% END