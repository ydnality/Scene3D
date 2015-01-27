%s_spatially_variant_wavefront_coeffs.m
%
%Estimate how each coeffs is affected by depth and field

% MP Vistasoft, 2014

%% Initialize Iset
s_initISET
wave = 400:50:600;

%Refractive indices
n_ob = 1 ; %object space
n_im = 1; %image space

%% 1.1 Create a multiple lens system declaring camera properties


% Declare film
% filmPosition = [0 0 51.2821	];  % Good for 2Elens
filmPosition = [0 0 37.4];        % Good for dgauss.50mm.  True focal about 37.3mm
% filmPosition = [0 0 107]; 
filmDiag = 3;  % Millimeters
filmSize = [filmDiag/sqrt(2) filmDiag/sqrt(2)];
resolution =  [300 300 length(wave)];
film = filmC('position', filmPosition, 'size', filmSize, 'wave', wave, 'resolution', resolution);
%
% lensFile = fullfile(s3dRootPath, 'data', 'lens', '2ElLens.mat');
lensFile = fullfile(s3dRootPath, 'data', 'lens', 'dgauss.50mm.mat');
import = load(lensFile,'lens');
lens = import.lens;
lens.apertureMiddleD = 4;
lens.set('wave', wave);

% Matrix of n for each surface element.
% Apertures are 0
n = lens.get('nArray');


%% 1.2 Declare point source in object space  and refractive indices of both object and image space

% The center of the first camera aperture is usually at [0 0 0].
% The object positions data are in -z.  We are using a right-handed
% coordinate system. 

% Specify the point source

% Angulare eccentricity (in degree)
angle_field=[0:2:30]; %�
% Distance 
% ps_dist=[ 1e3, 0.25*1e4, 0.5*1e4, 0.75*1e4, 1e4, 1e5, 1e6, 1e7, 1e8]; %distance
ps_dist=[ 1e2 1e3, 0.25*1e4, 0.5*1e4, 0.75*1e4, 1e4, 1e5]; %distance

for fi=1:length(angle_field)
    for di=1:length(ps_dist)
        % Create point source
        ps_x=ps_dist(di).*tan(angle_field(fi)/180*pi);
        pointSource = [ps_x 0 -ps_dist(di)];  % A very distant point.
        %Create camera
        psfCamera1 = psfCameraC('lens', lens, 'film', film, 'pointSource', pointSource);
        % Build Black Box Model
        psfCamera1.bbmCreate(n_ob,n_im);
        % Get coeffs
        defocusO=psfCamera1.get('bbm','defocus','#wave');        
        Delta_shift(fi,di,:) =psfCamera1.get('bbm','defocusshift');
        paO=psfCamera1.get('bbm','primaryaberration','#wave');
        % Store wavelength-dependent coeffs
        def(fi,di,:)=defocusO.W20;
        spher(fi,di,:)=paO.W40; %spherical aberration
        petzval(fi,di,:)=paO.W20; %petzval
        astig(fi,di,:)=paO.W22; %astigm
        coma(fi,di,:)=paO.W31; %coma
        dist(fi,di,:)=paO.W11; %distorsion
    end   
    
end


%% PLOT

%Specify a wavelength

wave0=550; %nm
indW0=find(wave==wave0);


% SURF
[pd,af]=meshgrid((ps_dist*1e-3),angle_field);
xL=['Field eccentricity [�]'];yL=['Distance  [m]'];zL=['Wavefront Coeff [#wave]'];
yLl=['Distance  log[m]'];
figure
surf(af,log(pd),def(:,:,indW0))
xlabel(xL),ylabel(yLl),zlabel(zL)
title(['Defocus: spatially dependence',' for lambda:',num2str(wave0),' nm'])
figure
surf(af,log(pd),spher(:,:,indW0))
xlabel(xL),ylabel(yLl),zlabel(zL)
title(['Spherical Aberration: spatially dependence',' for lambda:',num2str(wave0),' nm'])
figure
surf(af,(pd),petzval(:,:,indW0))
xlabel(xL),ylabel(yL),zlabel(zL)
title(['Petzval curvature: spatially dependence',' for lambda:',num2str(wave0),' nm'])
figure
surf(af,(pd),astig(:,:,indW0))
xlabel(xL),ylabel(yL),zlabel(zL)
title(['Astigmatism: spatially dependence',' for lambda:',num2str(wave0),' nm'])
figure
surf(af,(pd),coma(:,:,indW0))
xlabel(xL),ylabel(yL),zlabel(zL)
title(['Coma: spatially dependence',' for lambda:',num2str(wave0),' nm'])
figure
surf(af,(pd),dist(:,:,indW0))
xlabel(xL),ylabel(yL),zlabel(zL)
title(['Distorsion: spatially dependence',' for lambda:',num2str(wave0),' nm'])


% %IMAGESC
% pd_v=(ps_dist*1e-3); af_v=angle_field;
% 
% xL=['Field eccentricity [�]'];yL=['Distance  [m]'];zL=['Wavefront Coeff [#wave]'];
% yLl=['Distance  log[m]'];
% figure
% imagesc(af_v,log(pd_v),def(:,:,indW0))
% xlabel(xL),ylabel(yLl),zlabel(zL)
% title(['Defocus: spatially dependence',' for lambda:',num2str(wave0),' nm'])
% figure
% imagesc(af_v,log(pd_v),spher(:,:,indW0))
% xlabel(xL),ylabel(yLl),zlabel(zL)
% title(['Spherical Aberration: spatially dependence',' for lambda:',num2str(wave0),' nm'])
% figure
% imagesc(af_v,(pd_v),petzval(:,:,indW0))
% xlabel(xL),ylabel(yL),zlabel(zL)
% title(['Petzval curvature: spatially dependence',' for lambda:',num2str(wave0),' nm'])
% figure
% imagesc(af_v,(pd_v),astig(:,:,indW0))
% xlabel(xL),ylabel(yL),zlabel(zL)
% title(['Astigmatism: spatially dependence',' for lambda:',num2str(wave0),' nm'])
% figure
% imagesc(af_v,(pd_v),coma(:,:,indW0))
% xlabel(xL),ylabel(yL),zlabel(zL)
% title(['Coma: spatially dependence',' for lambda:',num2str(wave0),' nm'])
% figure
% imagesc(af_v,(pd_v),dist(:,:,indW0))
% xlabel(xL),ylabel(yL),zlabel(zL)
% title(['Distorsion: spatially dependence',' for lambda:',num2str(wave0),' nm'])




