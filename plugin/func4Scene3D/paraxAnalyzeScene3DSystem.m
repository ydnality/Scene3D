% analyze the  SCENE 3D system through paraxial optics
function [result]=paraxAnalyzeScene3DSystem(type,Syst,varargin)


%INPUT
%Syst:system of Scene 3D
%type: specify the type of object to analyze {'lens';'surfaceArray';'specify'}
%       in case of 'specify' to size (varargin,1)=2 were film and
%       poinSource has to be specify
%OUTPUT
%result: structure of output

%% CHECK INPUT

switch type
    case {'Lens';'lens'}
        %Get experimental conditions
        unit='mm'; %unit
        wave=Syst.wave*1e-6; % in mm
        nw=length(wave); %num wavelength
        
        nelem=length(Syst.surfaceArray);
        %Initialize some vector
        N=ones(nw,nelem); 
        %Useful parameter
        inD=1;
        %% Get the parameter to build the Optical System
        for ni=1:nelem
            %Get the structure
            S=Syst.surfaceArray(ni);
            if all(S.n==0)
                surftype{ni}='diaphragm';  
                if (S.apertureD==Syst.apertureMiddleD) %Check if the aperture size is changed
                    Diam(ni)=S.apertureD; %aperture diameter
                else
                    Diam(ni)=Syst.apertureMiddleD; %set aperture change
                end
                %save indices of the aperture
                indDiaph(inD)=ni;
                inD=inD+1;
                
                if ni>1
                    N(:,ni)=N(:,ni-1); %refractive indices
               end
            else
                surftype{ni}='refr';
                N(:,ni)=S.n';           %refractive indices                
                Diam(ni)=S.apertureD; %aperture diameter
            end
            Radius(ni)=S.sRadius; %radius of curvature
            posZ(ni)=S.get('zintercept');
        end
        % Set new origin as the First Surface
        PosZ=posZ-posZ(1);
        
        %% Build several surface
        for k=1:length(Radius)
            R(k)=Radius(k);
            z_pos(k)=PosZ(k);
            n(:,k)=N(:,k);
            diam(k)=Diam(k);
            switch surftype{k}
                case {'refr'}          
                    surf{k}=paraxCreateSurface(z_pos(k),diam(k),unit,wave,surftype{k},R(k),n(:,k));
                case {'diaphragm','diaph'}
                    surf{k}=paraxCreateSurface(z_pos(k),diam(k),unit,wave,surftype{k});
                case {'film'}
            end
        end


        %% CREATE OPTICAL SYSTEM
        if nargin>2
            n_ob=varagin{1};n_im=varargin{2};
        else            
            n_ob=1; n_im=1;
        end
        [OptSyst]=paraxCreateOptSyst(surf,n_ob,n_im,unit,wave);
        
        %% SET OUTPUT
        % parameter to coord change
        z0=OptSyst.cardPoints.lastVertex;
        
        result.focallength=OptSyst.cardPoints.fi; %focal lenght of the system
        % Cardinal Point
        result.cardinalPoint.ImageSpace.focalPoint=OptSyst.cardPoints.dFi; %Focal point in the image space
        result.cardinalPoint.ImageSpace.principalPoint=OptSyst.cardPoints.dHi; % Principal point in the image space
        result.cardinalPoint.ImageSpace.nodalPoint=OptSyst.cardPoints.dNi; % Nodal point in the image space
        result.cardinalPoint.ObjectSpace.focalPoint=OptSyst.cardPoints.dFo-z0; %Focal point in the object space
        result.cardinalPoint.ObjectSpace.principalPoint=OptSyst.cardPoints.dHo-z0; % Principal point in the object space
        result.cardinalPoint.ObjectSpace.nodalPoint=OptSyst.cardPoints.dNo-z0; % Nodal point in the object space
        result.abcdMatrix=OptSyst.matrix.abcd; % The 4 coefficients of the ABCD matrix of the overall system
        result.focalPlaneRadius=OptSyst.Petzval.radius; % radius of curvature of focal plane
    case {'surfaceArray';'surface Array';'psfCamera';' psf camera'} 
        %Get inputs
        lens=Syst.lens;
        film=Syst.film;
        pSource=Syst.pointSource;
        %Get experimental conditions
        unit='mm'; %unit
        wave=lens.wave*1e-6; % in mm
        nw=length(wave); %num wavelength
        
        nelem=length(lens.surfaceArray);
        %Initialize some vector
        N=ones(nw,nelem); 
        %Useful parameter
        inD=1;
        %% Get the parameter to build the Optical System
        for ni=1:nelem
            %Get the structure
            S=lens.surfaceArray(ni);
            if all(S.n==0)
                surftype{ni}='diaphragm';
                %save indices of the aperture
                indDiaph(inD)=ni;
                inD=inD+1;
                if (S.apertureD==Syst.lens.apertureMiddleD) %Check if the aperture size is changed
                    Diam(ni)=SapertureD; %aperture diameter
                else
                    Diam(ni)=Syst.lens.apertureMiddleD; %set aperture change
                end
                if ni>1
                    N(:,ni)=N(:,ni-1); %refractive indices
                end
            else
                surftype{ni}='refr';
                N(:,ni)=S.n';           %refractive indices
                Diam(ni)=S.apertureD; %aperture diameter
            end
            Radius(ni)=S.sRadius; %radius of curvature            
            posZ(ni)=S.get('zintercept');
        end
        % Set new origin as the First Surface
        PosZ=posZ-posZ(1);
        
        %% Build several surface
        for k=1:length(Radius)
            R(k)=Radius(k);
            z_pos(k)=PosZ(k);
            n(:,k)=N(:,k);
            diam(k)=Diam(k);
            switch surftype{k}
                case {'refr'}          
                    surf{k}=paraxCreateSurface(z_pos(k),diam(k),unit,wave,surftype{k},R(k),n(:,k));
                case {'diaphragm','diaph'}
                    surf{k}=paraxCreateSurface(z_pos(k),diam(k),unit,wave,surftype{k});
                case {'film'}
            end
        end


        %% CREATE OPTICAL SYSTEM
        if nargin>2
            n_ob=varagin{1};n_im=varargin{2};
        else            
            n_ob=1; n_im=1;
        end
        [OptSyst]=paraxCreateOptSyst(surf,n_ob,n_im,unit,wave);
        
        %% CREATE IMAGING SYSTEM
        %film object
        film_zpos=film.position(3)+OptSyst.cardPoints.lastVertex;
        profile='flat'; resFilm=film.resolution(1:2);pixel_pitch=film.size; %um x um
        [filmObj]=paraxCreateFilm(0,profile,resFilm,pixel_pitch,unit);
        %Imaging system
        [ImagSyst]=paraxCreateImagSyst(OptSyst,filmObj,film_zpos,[0 0]);
        % Point source object
        ps_height=sqrt(pSource(1).^2+pSource(2).^2); %distance of the point source to optical axis
        if not(pSource(2)==0) && not(pSource(1)==0)
            ps_angle=atan(pSource(2)/pSource(1)); %angle subtended by the point source and the x-axis in the object plane
        else
            if (pSource(2)==0)
                ps_angle=0;
            else
                ps_angle=pi/2;
            end
        end
        ps_zpos=pSource(3)+OptSyst.cardPoints.lastVertex; %point source position along the optical axis
        [pSourceObj]=paraxCreateObject(ps_zpos,ps_height,profile,unit);
        %Add to the Imaging System
        [ImagSyst]=paraxAddObject2ImagSyst(ImagSyst,pSourceObj);
        %% SET OUTPUT
        % parameter to coord change
        z0=OptSyst.cardPoints.lastVertex;
        
        result.focallength=ImagSyst.cardPoints.fi; %focal lenght of the system
        % Cardinal Point
        result.cardinalPoint.ImageSpace.focalPoint=ImagSyst.cardPoints.dFi; %Focal point in the image space
        result.cardinalPoint.ImageSpace.principalPoint=ImagSyst.cardPoints.dHi; % Principal point in the image space
        result.cardinalPoint.ImageSpace.nodalPoint=ImagSyst.cardPoints.dNi; % Nodal point in the image space
        result.cardinalPoint.ObjectSpace.focalPoint=ImagSyst.cardPoints.dFo-z0; %Focal point in the object space
        result.cardinalPoint.ObjectSpace.principalPoint=ImagSyst.cardPoints.dHo-z0; % Principal point in the object space
        result.cardinalPoint.ObjectSpace.nodalPoint=ImagSyst.cardPoints.dNo-z0; % Nodal point in the object space
        % Image formation
        result.abcdMatrix=ImagSyst.matrix.abcd; % The 4 coefficients of the ABCD matrix of the overall system
        result.focalPlaneRadius=ImagSyst.Petzval.radius; % radius of curvature of focal plane
        result.Fnumber=ImagSyst.object{end}.Radiance.Fnumber.eff; %Effective F number
        NA=n_im*sin(atan(ImagSyst.object{end}.Radiance.ExP.diam(:,1)./(ImagSyst.object{end}.ConjGauss.z_im-mean(ImagSyst.object{end}.Radiance.ExP.z_pos,2))));
        result.numericalAperture=NA; %effective numerical aperture
        magn_lateral=ImagSyst.object{end}.ConjGauss.m_lat; %
        magn_angular=ImagSyst.object{end}.ConjGauss.m_ang;
        result.magnification.lateral=magn_lateral; %lateral magnification
        result.magnification.angular=magn_angular; %lateral magnification
        % Image Point
        iPoint_zpos=ImagSyst.object{end}.ConjGauss.z_im-z0; %image point z position
        iPoint_heigh=ps_height.*magn_lateral; % image point distance from the optical axis
        iPoint_xpos=iPoint_heigh*cos(ps_angle); % image point x position
        iPoint_ypos=iPoint_heigh*sin(ps_angle); %image point y position
        result.imagePoint.position=[iPoint_xpos,iPoint_ypos,iPoint_zpos]; % x,y,z  coords for image point
        % Aperture and Pupils
        %Exit Pupil
        result.ExitPupil.zpos=mean(ImagSyst.object{end}.Radiance.ExP.z_pos,2);
        result.ExitPupil.diam=ImagSyst.object{end}.Radiance.ExP.diam(:,1)-ImagSyst.object{end}.Radiance.ExP.diam(:,2);
        %Entrance Pupil
        result.EntrancePupil.zpos=mean(ImagSyst.object{end}.Radiance.EnP.z_pos,2);
        result.EntrancePupil.diam=ImagSyst.object{end}.Radiance.EnP.diam(:,1)-ImagSyst.object{end}.Radiance.EnP.diam(:,2);
        % Wavefront Aberration Coefficients
        result.aberration.paCoeff=paUnit2NumWave(ImagSyst.object{end}.Wavefront.PeakCoeff);
    case {'specify'}
        % Get input
        lens=Syst;
        if nargin==4
            film=varargin{1};
            pSource=varargin{2};
        else
            error('Not enough inputs')
        end
        % Create psfCamera
        ppsfCamera = ppsfCameraC('lens', lens, 'film', film, 'pointSource', pSource);
        [result]=paraxAnalyzeScene3DSystem('psfCamera',ppsfCamera);
        
    otherwise
        error (['Not accepted ',type,' as system type']);

end
