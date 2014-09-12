function [ AInterp A1stInterp A2ndInterp ] = interpolateA( obj, wantedPSLocation, wavelength )
%INTERPOLATEA Summary of this function goes here
%   Detailed explanation goes here 
%   wavelength is ignored for now.  TODO: deal with wavelength and depth

% This is going to be the basis of the 'get' part of the VOLT class, when
% we return an interpolated linear transformation
AInterp = zeros(4,4);
A1stInterp = zeros(4,4);
A2ndInterp = zeros(4,4);

AComplete = obj.get('ACollection');
A1stComplete = obj.get('A1stCollection');
A2ndComplete = obj.get('A2ndCollection');

pSY = obj.get('fieldPositions');

for i = 1:4
    for j = 1:4
        coefValues = AComplete(i,j,:);
        coefValues = coefValues(:);
        yi = interp1(pSY,coefValues, wantedPSLocation);
        AInterp(i,j) = yi;
        
        coefValues = A1stComplete(i,j,:);
        coefValues = coefValues(:);
        yi = interp1(pSY,coefValues, wantedPSLocation);
        A1stInterp(i,j) = yi;
        
        coefValues = A2ndComplete(i,j,:);
        coefValues = coefValues(:);
        yi = interp1(pSY,coefValues, wantedPSLocation);
        A2ndInterp(i,j) = yi;
    end
end   

end

