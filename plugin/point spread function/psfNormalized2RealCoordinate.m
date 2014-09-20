%Conver the IMAGE COORDINATE from NORMALIZED to REAL 

%normalized coordinate (xim,yim) are linked to the image plane 'real'
% coordinate (Xi,Yi) from the following relation: 
% Xi=xp*lambda/NA
% Yi=yp*lambda/NA
%






function [Xi,Yi]=psfNormalized2RealCoordinate(X,Y,wave,NA)

%%INPUT
%X: x normalized coord in the image plane [1xNx]  Nx samples along X
%Y: y nomalized coord in the imageplane [1xNy]  Ny samples along y
%wave: wavelength sample [Mx1] M number of wavelength samples
%NA: Numerical aperture[scalar or Mx1]

%OUTPUT
%Xi=x coord in the image plane. [MxNx] wavelength dependent
%Yi= y coord in the image plane [MxNy] wavelength dependent

%NOTE
% all the distances have to be in the same unit e.g. [mm]
% 
% %% CHECK sampling
% if numel(delta_angle)==1
%     DangleX=delta_angle;
%     DangleY=delta_angle; 
%     
% elseif numel(delta_angle)==2
%     DangleX=delta_angle(1);
%     DangleY=delta_angle(2);
% else
%     error('Not accepted more that 2 value for sampling')
% end


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
%     Xi=X.*ratio/(2*pi);
%     Yi=Y.*ratio/(2*pi);
    Xi=X.*ratio;    
    Yi=Y.*ratio;
end

    