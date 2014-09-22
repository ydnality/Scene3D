


function [ImagSyst,vararging]=paraxOpt2Imag(OptSys,F,pSource,unit)


% Build imaging system composed by Optical System (OptSyst) + Film + point Source

% INPUT
% OptSys: struct of the optical system,  analisis through parax optics
%
% F: struct   .z= position
%                .res= resolution [Nx,Ny]
%                .pp= pixel pitch
% pSource:  [x y z ]
%
% unit: 'mm'

% OUTPUT
% ImagSyst: struct of the imaging system,  analisis through parax optics
%


%% Default parameters

%% CREATE IMAGING SYSTEM
%film object
film_zpos=F.z;
profile='flat'; resFilm=F.res; pixel_pitch=F.pp; %um x um


% Point source object
%         ps_height=sqrt(pSource(1).^2+pSource(2).^2); %distance of the point source to optical axis
%         if not(pSource(2)==0) && not(pSource(1)==0)
%             ps_angle=atan(pSource(2)/pSource(1)); %angle subtended by the point source and the x-axis in the object plane
%         else
%             if (pSource(2)==0)
%                 ps_angle=0;
%             else
%                 ps_angle=pi/2;
%             end
%         end
%         ps_zpos=pSource(3)+OptSyst.cardPoints.lastVertex; %point source position along the optical axis
%         

[filmObj]=paraxCreateFilm(0,profile,resFilm,pixel_pitch,unit);
%Imaging system
[ImagSyst]=paraxCreateImagSyst(OptSys,filmObj,film_zpos,[0 0]);
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
ps_zpos=pSource(3)+OptSys.cardPoints.lastVertex; %point source position along the optical axis
[pSourceObj]=paraxCreateObject(ps_zpos,ps_height,profile,unit);
%Add to the Imaging System
[ImagSyst]=paraxAddObject2ImagSyst(ImagSyst,pSourceObj);