

function [iPoint]=findImagePoint(obj,pSource,n_ob,n_im)


% Find the image point for a given object point thought the 
% imaging system specify by the principal points (in object space and image
% space), focal lenght

% INPUT
% pSource: point source coordinate [ x y z]
% n_ob: refractive index in the object space
% n_im: refractive index in the image space

%OUTPUT
%[iPoint]= [x y z]



%% GET POINT SOURCE POLAR COORDINATEs
[ps_height,ps_angle,ps_zpos]=coordCart2Polar3D(pSource(1),pSource(2),pSource(3)); %get image coordinate in polar coordinate

%% GET LENS PARAMETERs

bbm=obj.BBoxModel;
%principal point in the object space
H_obj=obj.bbmGetValue('objectprincipalpoint');
% Hobj=result2.cardinalPoint.ObjectSpace.principalPoint; 
%principal point in the image space
H_im=obj.bbmGetValue('imageprincipalpoint'); 
% Him=result2.cardinalPoint.ImageSpace.principalPoint; 

%Focal length
focalLength=obj.bbmGetValue('effectivefocallength'); 
% focalLength=result2.focallength; %effective focal length


%% IMAGE FORMATION
% Object-principal plane distance
d_ob=H_obj-ps_zpos; % distance between object point and  related principal point (Hob)

% equation    n_ob/d_ob  + n_im/d_im =1/efl
P=1./focalLength; %optical power
T1=P-n_ob./d_ob;
d_im=n_im./T1; % distance between image point and related principal point (H_im)

%Lateral magnification
%
%  m=-  (d_im/d_ob) (n_ob/n_im)
m_lat=-(d_im./d_ob).*(n_ob./n_im);

%% Image point
ip_zpos=H_im+d_im; %image point position along the optical axis

ip_height=ps_height.*m_lat; %image point distance from the optical axis
% the angle is unchanged
% ip_xpos=ip_height.*cos(ps_angle); % x coordinate
% ip_ypos=ip_height.*sin(ps_angle); % x coordinate

%% SET OUTPUT
[iPoint(:,1),iPoint(:,2),iPoint(:,3)]=coordPolar2Cart3D(ip_height,ps_angle,ip_zpos);
% iPoint=[ip_xpos ip_ypos ip_zpos]; 