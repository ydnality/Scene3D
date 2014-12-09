function [ reflectanceEst, meas2BasisCoeffs, coeffs, predResp ] = linearFindReflectanceV2( measVals, cameraResp, cameraGain, cameraOffset, illuminant, basisFcns, lambda )

% [ reflectanceEst, meas2BasisCoeffs, coeffs, predResp ] = linearFindReflectanceV2( measVals, cameraResp, cameraGain, cameraOffset, illuminant, basisFcns, lambda )
%
% This function computes the surface reflectance estimates given camera and
% illumination parameters. It is assumed that the camera has "nFilter"
% different spectral channels and the scene is observed under "nChannel"
% illuminants. 
%
% This function estimates the reflectance by approximating the
% spectral reflectance with a small number of basis functions and computing
% a least-squares fit of these coefficients to acquited data.
%
% Inputs:
%   measVals - a nFilter x nChannel x nSamples matrix with the pixel values
%   cameraResp - a nWaves x nFilter matrix with the spectral responsivity
%   of each camera channel
%   cameraGain - a nFilter x nChannel matrix with gains applied to each of
%   the channel-illuminant combinations
%   cameraOffset - a nFilter x nChannel matrix with offsets for each
%   channel-illuminant combination
%   illuminant - a nWaves x nChannel matrix with illumination intensity (in
%   photons)
%   basisFcns - a nWaves x nBasis matrix with basis function
%   representations for the spectral reflectance.
%   lambda - a scalar describing the smoothness parameter
%
% Outputs:
%   reflectanceEst - a nWaves x nSamples matrix with spectral reflectance
%   estimates
%   meas2BasisCoeffs - a nBasis x (nFilters*nChannels) matrix that converts
%   the measured values to basis function coefficients
%   coeffs - a nBasis x nSamples matrix with basis function weights for
%   each reflectance estimate
%   predResp - a nFilter x nChannel x nSamples matrix with modeled pixel
%   responses.
%
% Copyright, Henryk Blasinski 2014

nBasis = size(basisFcns,2);
nChannels = size(illuminant,2);
nFilters = size(cameraResp,2);
nSamples = size(measVals,3);
nWaves = size(illuminant,1);

% Correct the measured values for the camera offset
corrMeasVals = measVals - repmat(cameraOffset,[1,1,nSamples]);
b = reshape(corrMeasVals,nChannels*nFilters,nSamples);

% Basis vectors
basisVectors = zeros(nChannels*nFilters,nBasis);
for i=1:nBasis
    tmp = cameraGain.*(cameraResp'*diag(basisFcns(:,i))*illuminant);
    basisVectors(:,i) = tmp(:);
end

% Create a penalty on roughness
Rm = [diag(ones(nWaves-1,1)) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) diag(ones(nWaves-1,1))];
Rm = Rm*basisFcns;

% Create the image formation model matrix
A = [basisVectors; sqrt(lambda)*Rm];

% Compute the inversion
meas2BasisCoeffs = pinv(A);
meas2BasisCoeffs = meas2BasisCoeffs(:,1:nChannels*nFilters);

coeffs = meas2BasisCoeffs*b;
reflectanceEst = basisFcns*coeffs;

predResp = basisVectors*coeffs;
predResp = reshape(predResp,[nFilters, nChannels, nSamples]);
predResp = predResp + repmat(cameraOffset,[1, 1, nSamples]);
    
end

