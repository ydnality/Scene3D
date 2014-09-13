function [ AInterp A1stInterp A2ndInterp ] = interpolateA( obj, wantedPSLocation, wavelength )
%INTERPOLATEA Summary of this function goes here
%   Detailed explanation goes here 
%   wavelength is ignored for now.  TODO: deal with wavelength and depth
%
%   For now - we will assume that x coordinates are 0 and interpolation can
%   only be done on the y axis and depth.  TODO: allow for rotations for
%   anything out of this plane.

% This is going to be the basis of the 'get' part of the VOLT class, when
% we return an interpolated linear transformation
AInterp = zeros(4,4);
A1stInterp = zeros(4,4);
A2ndInterp = zeros(4,4);

AComplete = obj.get('ACollection');
A1stComplete = obj.get('A1stCollection');
A2ndComplete = obj.get('A2ndCollection');

numDepths = obj.get('numDepths');
numPositions = obj.get('numFieldPositions');
pSY = obj.get('fieldPositions');
pSZ = obj.get('depths');

[meshY, meshZ] = meshgrid(pSY,pSZ);


for i = 1:4
    for j = 1:4
        coefValues = AComplete(i,j,:,:);
        %coefValues = coefValues(:);
        %coefValues dimensions: (1,1, #fieldPositions, #depths);
        coefValues = reshape(coefValues, numDepths, numPositions);
        yi = interp2(meshY, meshZ, coefValues, wantedPSLocation(2), wantedPSLocation(3));
        AInterp(i,j) = yi;
        
        coefValues = A1stComplete(i,j,:,:);
        %coefValues = coefValues(:);
        coefValues = reshape(coefValues, numDepths, numPositions);
        yi = interp2(meshY, meshZ, coefValues, wantedPSLocation(2), wantedPSLocation(3));
        A1stInterp(i,j) = yi;
        
        coefValues = A2ndComplete(i,j,:,:);
        %coefValues = coefValues(:);
        %yi = interp1(pSY,coefValues, wantedPSLocation);
        coefValues = reshape(coefValues, numDepths, numPositions);
        yi = interp2(meshY, meshZ, coefValues, wantedPSLocation(2), wantedPSLocation(3));
        A2ndInterp(i,j) = yi;
    end
end   

end

