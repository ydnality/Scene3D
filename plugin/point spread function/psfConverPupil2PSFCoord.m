%Conver the normalized coordinate of the PUPIL to the PSF in the image plane
% 
% Pupil normalized coordinate (xp,yp) are linked to the image plane 'real'
% coordinate (Xi,Yi) from the following relation: 
% Xi=xp*lambda/NA
% Yi=yp*lambda/NA
%




% In pupil plane (Xp,Yp)->In image plane (Xi,Yi)  Xi=Xp *pi/(lambda*
% effectiveF-number), same for Yi
%Effective F number is the ratio between the ExP-GaussianImagePoint and the
%ExP Diameter


function [Xi,Yi]=psfConverPupil2PSFCoord(Xp,Yp,delta_angle,wave,NA)

%%INPUT
%Xp: x normalized coord in the pupil plane [1xNx]  Nx samples along X
%Yp: y nomalized coord in the pupil plane [1xNy]  Ny samples along y
%delta_angle: sampling of the fourier transform, scalar :-> same sampling
%for X and Y, [1x2]:-> (1) sampling for x, (2) sampling for y
%wave: wavelength sample [Mx1] M number of wavelength samples
%NA: Numerical aperture[scalar or Mx1]

%OUTPUT
%Xi=x coord in the image plane. [MxNx] wavelength dependent
%Yi= y coord in the image plane [MxNy] wavelength dependent

%NOTE
% all the distances have to be in the same unit e.g. [mm]

%% CHECK sampling
if numel(delta_angle)==1
    DangleX=delta_angle;
    DangleY=delta_angle; 
    
elseif numel(delta_angle)==2
    DangleX=delta_angle(1);
    DangleY=delta_angle(2);
else
    error('Not accepted more that 2 value for sampling')
end


%% EARLY PARAMETERs

nw=size(wave,1);


%% CHECK if effective F number is wavelength dependent



for li=1:nw
    if numel(NA)==1
        ratio=(wave(li,1)./NA); %scaling factor
    elseif numel(NA)==nw
        ratio=(wave(li,1)./NA(li,1)); %scaling factor
    else
        error('effective F number and wavelength DO NOT MATCH !')
    end    
    %Rescale the variable
    Xi=Xp.*ratio.*DangleX;
    Yi=Yp.*ratio.*DangleY;
end

    