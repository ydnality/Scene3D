% Find the image point for a given object point thought the 
% imaging system specify by the principal points (in object space and image
% space), focal lenght

function [iPoint]=findImagePoint(pSource,n_ob,H_ob,H_im,n_im,focalLength)

% INPUT
% pSource: point source coordinate [ x y z]
% n_ob: refractive index in the object space
% H_ob: principal point position along the optical axis for the object % space
% H_im: principal point position along the optical axis for the image % space
% n_im: refractive index in the image space
% focalLength: focal lenght of the system

%OUTPUT
%[iPoint]= [x y z]



%% Object space
ps_zpos=pSource(3); %position  along the optical axis
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

d_ob=H_ob-ps_zpos; % distance between object point and  related principal point (Hob)

%% IMAGE FORMATION
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
ip_xpos=ip_height.*cos(ps_angle); % x coordinate
ip_ypos=ip_height.*sin(ps_angle); % x coordinate

%% SET OUTPUT
iPoint=[ip_xpos ip_ypos ip_zpos]; 