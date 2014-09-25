function [result,varargout] = paraxAnalyzeScene3DSystem(type,Syst,varargin)
% Analyze the  SCENE 3D system through paraxial optics
%
%  [result] = paraxAnalyzeScene3DSystem(type,Syst,varargin)
%
%
%INPUT
% type: specify the type of object to analyze {'lens';'surfaceArray';'specify'}
%       in case of 'specify' the to size (varargin,1)=2 were film and
%       poinSource has to be specify
% Syst:system of Scene 3D
%
%OUTPUT
%
%   result: A structure containing a variety of fields.  THe returned
%     fields can depend on which type of structure is sent in to the
%     routine. 
%
%   lens:
%
%   surface array:
%
%
%   all:  (lens, film, point source)
%
%
% MP Vistasoft 2014
  

%% TODO
%
%  Maybe call this bbmCreate(object, varargin)
%    type = object.type;
%

%% CHECK INPUT
% type = obj.type;
% 

type = ieParamFormat(type);

unit='mm'; %unit


%%        
switch type
    case 'lens'
        lens0=Syst;
        OptSys = paraxAnalyze(lens0);
        % Equivalent Black Box Model of lens
        lens0.set('black box model',OptSys);    
        % Set OUTPUT
        result=lens0;
       
    case {'psfcamera'} 

        % 
        % psfCamera = imgsysCompute(psfCamera);
        
        
        %Get inputs
        lens0=Syst.lens;
        film0=Syst.film;
        pSource=Syst.pointSource;
        
        % Equivalent Black Box Model of lens
        OptSys = paraxAnalyze(lens0);
        lens0.set('black box model',OptSys);    % Equivalent Black Box Model of lens
        
        % Build Imaging System composed by Optical System+Film+ pSource
        F.z=film0.position(3)+OptSys.cardPoints.lastVertex;
        F.res=film0.resolution(1:2);F.pp=film0.size; %um x um
        
        [ImagSyst]=paraxOpt2Imag(OptSys,F,pSource,unit); 
        
        [ps_heigth,ps_angle,ps_zpos]=coordCart2Polar3D(pSource(1),pSource(2),pSource(3)); %get image coordinate in polar coordinate
        
        pSpolar(1)=ps_heigth;pSpolar(2)=ps_angle;pSpolar(3)=ps_zpos;
        % Equivalent Black Box Model of Imaging System
        Syst.set('black box model',ImagSyst,pSpolar);
        
        % SET OUTPUT/s
        result=Syst;              
        % For more output
        if nargout>1
            varargout{1}=lens0;
        end
        
        
        
    case {'all'}
        % Suggest deleting this and making people use either lens or
        % psfCamera.  
        %
        % lens, film and point source are all provided
        
        % Get input
        lens=Syst;
        if nargin==4
            film=varargin{1};
            pSource=varargin{2};
        else
            error('Not enough inputs')
        end
        % Create psfCamera
        psfCamera = psfCameraC('lens', lens, 'film', film, 'pointSource', pSource);
        
        [result]=paraxAnalyzeScene3DSystem('psfCamera',psfCamera);
        % For more output
        if nargout>1
            [result,lens0]=paraxAnalyzeScene3DSystem('psfCamera',psfCamera);
            varargout{1}=lens0;
        else
            [result]=paraxAnalyzeScene3DSystem('psfCamera',psfCamera);
        end
        
    otherwise
        error (['Not accepted [',type,'] as system type']);

end
