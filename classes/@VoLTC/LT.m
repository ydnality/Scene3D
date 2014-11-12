function LT = LT(obj, pSource)
% Return a linear transform object for a given point from VoLT object
%
%  LT = VoLT.LT(pSource)
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

% Either interpolate or retrieve the matrices for the pSource location and
% all the wavelengths 
[AInterp, A1stInterp, A2ndInterp ] = obj.interpolateA(pSource, wave);

% Return the linear transform object for this point
LT = LTC('wave', wave, ...
    'AInterp', AInterp, ...
    'A1stInterp', A1stInterp, ...
    'A2ndInterp', A2ndInterp); 

end
