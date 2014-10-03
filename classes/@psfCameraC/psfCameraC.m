classdef psfCameraC <  handle
    % Create a point spread camera object
    %
    % Spatial units throughout are mm
    %
    % AL Vistasoft Copyright 2014
    % See Also: ppsfCameraC.
    
    % Figure out the relationship between these rays and the ppsfRays in
    % the ppsfCameraC.
    properties
        type = 'psfcamera';
        name = 'default camera';
        
        lens;
        film;
        pointSource;
        rays;
        BBoxModel;
        fftPSF;
    end
    
    methods (Access = public)
        
        %default constructor
        function obj = psfCameraC(varargin)
            % psfCameraC('lens',lens,'film',film,'point source',point);
            
            for ii=1:2:length(varargin)
                p = ieParamFormat(varargin{ii});
                switch p
                    case 'lens'
                        obj.lens = varargin{ii+1};
                    case 'film'
                        obj.film = varargin{ii+1};
                    case 'pointsource'
                        obj.pointSource = varargin{ii+1};
                    case {'blackboxmodel';'blackbox';'bbm'}
                       obj.BBoxModel = varargin{ii+1};
                    case {'fftpsf';'fftPSF'}
                        obj.lens = varargin{ii+1};
                    otherwise
                        error('Unknown parameter %s\n',varargin{ii});
                end
            end
            
        end
        
        function [val] = get(obj,param,varargin)
            % psfCamera.get('parameter name')
            % Start to set up the gets for this object
            val = [];
            param = ieParamFormat(param);
            switch param
                case {'wave','wavelength'}
                    val1 = obj.lens.wave(:);
                    val2 = obj.film.wave(:);
                    if isequal(val1,val2), val = val1; return;
                    else warning('Lens and film wavelength differ.  Using lens.');
                    end
                    val = val1;
                    
                case 'spacing'
                    % Millimeters per sample
                    r = obj.film.resolution(1);
                    s = obj.film.size(1);
                    val = s/r;
                case 'imagecentroid'
                    % obj.get('image centroid')
                    % x,y positions (0,0) is center of the image centroid.
                    % Used for calculating centroid of the psf
                    % Could use obj.film.image for the data, rather than oi
                    % p_renderOiMatlabToolFull
                    % Figure out center pos by calculating the centroid of illuminance image
                    flm = obj.film;
                    img = flm.image;  img = sum(img,3);

                    % Force to unit area and flip up/down for a point spread
                    img = img./sum(img(:));
                    img = flipud(img);
                    % vcNewGraphWin; mesh(img);

                    % Calculate the weighted centroid/center-of-mass
                    xSample = linspace(-flm.size(1)/2, flm.size(1)/2, flm.resolution(1));
                    ySample = linspace(-flm.size(2)/2, flm.size(2)/2, flm.resolution(2));
                    [filmDistanceX, filmDistanceY] = meshgrid(xSample,ySample);
                    
                    % distanceMatrix = sqrt(filmDistanceX.^2 + filmDistanceY.^2);
                    val.X = sum(sum(img .* filmDistanceX));
                    val.Y = sum(sum(img .* filmDistanceY));
                    
                 case {'blackboxmodel';'blackbox';'bbm'} % equivalent BLACK BOX MODEL
                     if nargin>2
                         fileType = varargin{1};  %witch field of the black box to get
                     else
                         error(['Specify also the field of the Black Box Model!'])
                     end
                     
                     if nargin>3                         
                        [val]=obj.bbmGetValue(fileType,varargin{2});
                     else                         
                        [val]=obj.bbmGetValue(fileType);
                     end
                     
%                     val=obj.bbmGetValue(obj.BBoxModel,fileType);
                        
                   case {'opticalsystem'; 'optsyst';'opticalsyst';'optical system structure'} 
                    % Get the equivalent optical system structure generated
                    % by Michael's script      
                    % Can be specify refractive indices for object and
                    % image space as varargin {1} and {2}
                    lens=obj.lens;
                    if nargin >2
                        n_ob = varargin{1};    n_im = varargin{2};
                        OptSyst = lens.bbmComputeOptSyst(n_ob,n_im);
                    else
                        OptSyst = lens.bbmComputeOptSyst();                    
                    end
                    val = OptSyst;
                    
                 case {'imagingsystem'; 'imgsyst';'imagingsyst';'imaging system structure'} 
                    % Get the equivalent imaging system structure generated
                    % by Michael's script      
                    % Can be specify refractive indices for object and
                    % image space as varargin {1} and {2}
                    %% Get inputs
                    lens=obj.lens;
                    film=obj.film;
                    pSource=obj.pointSource;
                    %COMPUTE OPTICAL SYSTEM
                    if nargin >2
                        n_ob = varargin{1};    n_im = varargin{2};
                        OptSyst = lens.get('optical system',n_ob,n_im);
                    else
                        OptSyst = lens.get('optical system');                  
                    end
                    unit=paraxGet(OptSyst,'unit');
                    % GET USEFUL PARAMETERs
                    lV=paraxGet(OptSyst,'lastvertex'); % last vertex of the optical system
%                     lV=0;                    
                    F.z=film.position(3)+lV;
                    F.res=film.resolution(1:2);F.pp=film.size; %um x um
                    
                    %CREATE an Imaging System
                   [ImagSyst]=paraxOpt2Imag(OptSyst,F,pSource,unit);       
                   
                   % SET OUTPUT
                    val = ImagSyst;
                    
                 case {'film'} 
                    % get the film structure   
                    val = obj.film;                    
                case {'pointsource';'psource';'lightsource'} 
                    % get the point source in the psfCamera   
                    val = obj.pointSource;
                case {'lens'} 
                    % get the lens the psfCamera   
                     val = obj.lens;
                case {'fftpsf';'psffft'} 
                    % get the fftPSF 
                     val = obj.fftPSF;   
                 case {'fftpsfmodulus';'psffftvalue'} 
                      % get the fftPSF modolus, 
                      % Specifying the wavelength if you want a specific PSF
                    if nargin>2
                        wave0=varargin{1};
                        waveV=obj.get('wave');
                        indW0=find(wave==wave0);
                        if isempty(indW0)
                            val=obj.fftPSF.abs;
                            warning (['The specified wavelength does match with the available ones: ',num2str(wave) ,' nm'])
                        else
                            val=obj.fftPSF.abs(:,:,indW0);
                        end
                    end
                    val=obj.fftPSF.abs;
                    case {'fftpsfcoordinate';'psffftcoord'} 
                      % get the fftPSF coord, 
                      % Specifying the wavelength if you want a specific PSF
                    if nargin>2
                        wave0=varargin{1};
                        waveV=obj.get('wave');
                        indW0=find(wave==wave0);
                        if isempty(indW0)
                            val.x=obj.fftPSF.x;
                            val.y=obj.fftPSF.y;
                            warning (['The specified wavelength does match with the available ones: ',num2str(wave) ,' nm'])
                        else
                            val.x=obj.fftPSF.x(:,:,indW0);
                            val.y=obj.fftPSF.y(:,:,indW0);
                        end
                    end
                    val.x=obj.fftPSF.x;
                    val.y=obj.fftPSF.y;
                otherwise
                    error('unknown parameter %s\n',param)
            end
            
        end
        
         function val = set(obj,param,val,varargin)
            % psfCamera.set('parameter name',value)
            % Start to set up the gets for this object
%             val = [];
            param = ieParamFormat(param);
            switch param
                case {'pointsource';'psource';'lightsource'} 
                    % set the point source in the psfCamera   
                    obj.pointSource= val;
                case {'lens'} 
                    % set the filmin the psfCamera   
                    obj.lens= val;
                case {'film'} 
                    % set the film in the psfCamera   
                    obj.film = val;
                case {'blackboxmodel';'blackbox';'bbm'};
                    %Get the parameters from the imaging system structure to build an  equivalent Black Box Model of the lens.
                    % The ImagSyst structure has to be built with the function 'paraxCreateImagSyst'
                    % Get 'new' origin for optical axis 
                    % INPUT
                    % val= ImagSyst struct
                    % varargin {1}: polar coordinate of pointSource [ro, theta, z]
                    %
                    % MP Vistasoft 2014
                    %
                    ImagSyst=val;
                    psPolar=varargin{1};
                    z0 = paraxGet(ImagSyst,'lastvertex');
                    % Variable to append
                    efl=paraxGet(ImagSyst,'effectivefocallength'); %focal lenght of the system
                    obj=obj.bbmSetField('effectivefocallength',efl);
                    pRad = paraxGet(ImagSyst,'effectivefocallength'); % radius of curvature of focal plane
                    obj=obj.bbmSetField('focalradius',pRad);
                    Fi=paraxGet(ImagSyst,'imagefocalpoint')-z0;     %Focal point in the image space
                    obj=obj.bbmSetField('imagefocalpoint',Fi);
                    Hi=paraxGet(ImagSyst,'imageprincipalpoint')-z0; % Principal point in the image space
                    obj=obj.bbmSetField('imageprincipalpoint',Hi);
                    Ni=paraxGet(ImagSyst,'imagenodalpoint')-z0;     % Nodal point in the image space
                    obj=obj.bbmSetField('imagenodalpoint',Ni);
                    Fo=paraxGet(ImagSyst,'objectfocalpoint')-z0; %Focal point in the object space
                    obj=obj.bbmSetField('objectfocalpoint',Fo);
                    Ho=paraxGet(ImagSyst,'objectprincipalpoint')-z0; % Principal point in the object space
                    obj=obj.bbmSetField('objectprincipalpoint',Ho);
                    No=paraxGet(ImagSyst,'objectnodalpoint')-z0; % Nodal point in the object space
                    obj=obj.bbmSetField('objectnodalpoint',No);
                    % abcd Matrix (Paraxial)
                    M = ImagSyst.matrix.abcd; % The 4 coefficients of the ABCD matrix of the overall system
                    obj=obj.bbmSetField('abcd',M);
                    
                    % IMAGE FORMATION                    
                    % Effective F number
                    Fnum=ImagSyst.object{end}.Radiance.Fnumber.eff; %effective F number
                    obj=obj.bbmSetField('fnumber',Fnum);
                    % Numerical Aperture
                    NA=ImagSyst.n_im.*sin(atan(ImagSyst.object{end}.Radiance.ExP.diam(:,1)./(ImagSyst.object{end}.ConjGauss.z_im-mean(ImagSyst.object{end}.Radiance.ExP.z_pos,2))));
                    obj=obj.bbmSetField('numericalaperture',NA);
                    %Field of View
                    FoV=ImagSyst.object{end}.Radiance.FoV;
                    obj=obj.bbmSetField('fieldofview',FoV);
                    % Lateral magnification
                    magn_lateral=ImagSyst.object{end}.ConjGauss.m_lat; %
                    obj=obj.bbmSetField('lateralmagnification',magn_lateral);                                   
                    % Exit Pupil
                    ExitPupil.zpos=mean(ImagSyst.object{end}.Radiance.ExP.z_pos,2)-z0;
                    ExitPupil.diam=ImagSyst.object{end}.Radiance.ExP.diam(:,1)-ImagSyst.object{end}.Radiance.ExP.diam(:,2);
                    obj=obj.bbmSetField('exitpupil',ExitPupil);
                    % Entrance Pupil
                    EntrancePupil.zpos=mean(ImagSyst.object{end}.Radiance.EnP.z_pos,2)-z0;
                    EntrancePupil.diam=ImagSyst.object{end}.Radiance.EnP.diam(:,1)-ImagSyst.object{end}.Radiance.EnP.diam(:,2);
                    obj=obj.bbmSetField('entrancepupil',EntrancePupil);                    
                    % Gaussian Image Point
                    iP_zpos=ImagSyst.object{end}.ConjGauss.z_im-z0; %image point z position
                    iP_h=psPolar(1).*magn_lateral;% image point distance from the optical axis
                    [iP(:,1),iP(:,2),iP(:,3)]=coordPolar2Cart3D(iP_h,psPolar(2),iP_zpos);
                    obj=obj.bbmSetField('gaussianimagepoint',iP);                    
                    % Aberration
                    % Primary Aberration
                    paCoeff=ImagSyst.object{end}.Wavefront.PeakCoeff;
                    obj=obj.bbmSetField('primaryaberration',paCoeff);
                    % Defocus
                    [obj_x,obj_y,obj_z]=coordPolar2Cart3D(psPolar(1),psPolar(2),psPolar(3)); 
%                     Obj.z=obj_z; Obj.y=obj_y;                    
                    Obj.z=obj_z+paraxGet(ImagSyst,'lastVertex'); 
                    Obj.y=sqrt(obj_x.^+obj_y.^2); % eccentricity (height)
                    [defCoeff] = paEstimateDefocus(ImagSyst,Obj,'best');
                    obj=obj.bbmSetField('defocus',defCoeff);
                                        
                    % REFRACTIVE INDEX
                    % object space
                    n_ob=ImagSyst.n_ob;
                    obj=obj.bbmSetField('n_ob',n_ob);
                    % image space
                    n_im=ImagSyst.n_im;
                    obj=obj.bbmSetField('n_im',n_im);
                    
                case {'fftpsf';'psffft'} 
                % get the fftPSF 
                  obj.fftPSF=val;   
                 case {'fftpsfmodulus';'psffftvalue'} 
                      % get the fftPSF modolus, 
                      % Specifying the wavelength if you want a specific PSF
                    if nargin>3
                        wave0=varargin{1};
                        waveV=obj.get('wave');
                        indW0=find(wave==wave0);
                        if isempty(indW0)
                            obj.fftPSF.abs=val;
                            warning (['The specified wavelength does match with the available ones: ',num2str(wave) ,' nm'])
                        else
                            obj.fftPSF.abs(:,:,indW0)=val;
                        end
                    end
                    obj.fftPSF.abs=val;
                    case {'fftpsfcoordinate';'psffftcoord'} 
                      % get the fftPSF coord, 
                      % Specifying the wavelength if you want a specific PSF
                    if nargin>3
                        wave0=varargin{1};
                        waveV=obj.get('wave');
                        indW0=find(wave==wave0);
                        if isempty(indW0)
                            obj.fftPSF.x=val.x;
                            obj.fftPSF.y=val.y;
                            warning (['The specified wavelength does match with the available ones: ',num2str(wave) ,' nm'])
                        else
                            obj.fftPSF.x(:,:,indW0)=val.x;
                            obj.fftPSF.y(:,:,indW0)=val.y;
                        end
                    end
                    obj.fftPSF.x=val.x;
                    obj.fftPSF.y=val.y;
                    
                otherwise
                    error('unknown parameter %s\n',param)
            end
            
        end
        
    end
    
end