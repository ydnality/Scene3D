function [Xi,Yi]=psfNormalized2RealCoordinate(X,Y,wave,NA)
%Convert from Normalized coordinates to real coordinates 
%
%  [Xi,Yi]=psfNormalized2RealCoordinate(X,Y,wave,NA)
%
% normalized coordinates in the image plane, (xim,yim) are linked to the
% image plane 'real' coordinate (Xi,Yi) from the following relation:
% 
%    Xi = xim * lambda/NA
%    Yi = yim * lambda/NA
%
%%INPUT
%   X: x normalized coord in the image plane [1xNx]  Nx samples along X
%   Y: y nomalized coord in the imageplane [1xNy]  Ny samples along y
%   wave: wavelength sample [Mx1] M number of wavelength samples
%   NA: Numerical aperture[scalar or Mx1]
%
%OUTPUT
%   Xi= x coord in the image plane. [MxNx] wavelength dependent
%   Yi= y coord in the image plane [MxNy] wavelength dependent
%
%NOTE
% The distance units should match, say they are all in mm
% 
% Example:
%
%
% See also:  psfEstimateImageCoord, psfPupil2PSF, ComputeFFTpsf.m
%
% MP Vistasoft Team, Copyright 2014

%% EARLY PARAMETERs

nw =size(wave,1);


%% CHECK if effective F number is wavelength dependent

for li=1:nw
    if numel(NA)==1
        ratio=(wave(li,1)./NA); %scaling factor
    elseif numel(NA)==nw
        ratio=(wave(li,1)./NA(li,1)); %scaling factor
    else
        error('Number of effective F numbers and number of wavelengths DO NOT MATCH !')
    end    
    Xi=X.*ratio;    
    Yi=Y.*ratio;
end

end
