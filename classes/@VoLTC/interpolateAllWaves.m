function [AInterp A1stInterp A2ndInterp ] = interpolateAllWaves(obj, wantedPSLocation )
%[AInterp A1stInterp A2ndInterp ] = interpolateAllWaves(obj, wantedPSLocation )
%   interpolates A matrices for all wavelengths.
%
%Inputs:
%   wantedPSLocation: a 1x3 vector containing the depth that the user wishes
%   to interpolate a set of linear matrices at
%
%Outputs:
%   AInterp: a 4 x 4 x w matrix.  Each 4 x 4 layer contains the 4 x 4
%   linear transform that summarizes light field transforms at that
%   particular point.  (will become 4 x 4 x n x w)
%
%   A1stInterp: The lens will be split into 2 layers: the first half
%   leading to the aperture.  This matrix will be a 4 x 4 x n x w matrix that
%   will be a collection of A matrices that summarizes the FIRST HALF of the
%   lens only.
%
%   A2ndInterp: a 4 x 4 x w matrix that contains a collection of A
%   matrices that summarizes the 2nd half of the lens.
%
% AL VISTASOFT 2014
%
% See Also: interpolateA

wave = obj.get('wave');
%full model from front-most element to back-most element
AInterp = zeros(4,4, length(wave));  
%half the model from the front-most element to the middle aperture
A1stInterp = zeros(4,4, length(wave));  
%the second half of the model from the middle aperture to the back-most element
A2ndInterp = zeros(4,4, length(wave));  

for w = 1:length(wave)
    [AInterp1Wave, A1stInterpCurrentWave, A2ndInterpCurrentWave ] = obj.interpolateA(wantedPSLocation, wave(w));
    AInterp(:,:,w) = AInterp1Wave;
    A1stInterp(:,:,w) = A1stInterpCurrentWave;
    A2ndInterp(:,:,w) = A2ndInterpCurrentWave;
end



end

