function [varargout] = bbmCreate(obj,varargin)
% Create the black box model of the lens 
%
% Various problems with the logic here, and we should re-write before long
% (BW).  See below.
%
% The black box model is a summary of the optical system that we use
% structure of SCENE3D to the Optical System structure to analyze data with
% paraxial approximation (first order optics) 
%
% This function does some strange things that we should think about.
%
% In this format, it creates the BBM with n_ob and n_im = 1
%  BBM = lens.bbmCreate()
%    or
%  BBM = lens.bbmCreate (n_ob,n_im)
%
% In this form, the lens call creates a psfCamera.  
%
%   psfCamera = lens.bbmCreate ('all',pointSource,film)
%       or
%   psfCamera = lens.bbmCreate ('all',pointSource,film,n_ob,n_im)
%
% The latter two calls seem weird to me (BW) and I think it should go away.
% I think we should only create the psfCamera explicitly and create the
% BBM.  This means removing the 'all' condition.
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
            psfFlag = 1; %YES psfCamera object to compute
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

% In this version MP sometimes returns a psfCamera and sometimes a black
% box.  Another reason to make this 'all' option go away.
if psfFlag
    psfCamera0 = psfCameraC('lens', lens0, 'film', film0, 'pointSource', pSource);
    psfCamera0.bbmCreate();
    varargout{1} = psfCamera0;
else
    varargout{1} =obj.get('bbm','all');
end

end

