function [c,oi] = centroid(obj,oi)
% Compute centroid of the oi illuminance distribution
%
%   [c,oi] = psfCamera.centroid(obj,oi)
%
% The coordinates are arranged with the middle of the film object as 0,0,
% and the number of samples are given by the film.size and film.resolution.
%
% The centroid is the intensity times the position.
% The position is the distance from the center (0,0).
%
% Inputs
%  obj: psfCamera
%  oi:  Precomputed optical image
%
% Returns
%  c.X and c.Y are the centroid coordinates
%  oi:  The optical image
%
% BW, Vistasoft, Copyright 2015

%% Perhaps we want to pass in the oi some time, for efficiency
if ~exist('oi','var')
    % Estimate the PSF and get the oi
    obj.estimatePSF;
    oi = obj.oiCreate(); 
end

%% Figure out center pos by calculating the centroid of illuminance image
img = oiGet(oi,'illuminance');

% Force to unit area and flip up/down for a point spread
img = img./sum(img(:));
img = flipud(img);    % Not sure why this is here
% vcNewGraphWin; mesh(img);

%% Calculate the weighted centroid/center-of-mass

% The units of the obj.film spacing are millimeters, I think (BW)
xSample = linspace(-obj.film.size(1)/2, obj.film.size(1)/2, obj.film.resolution(1));
ySample = linspace(-obj.film.size(2)/2, obj.film.size(2)/2, obj.film.resolution(2));
[DistanceX, DistanceY] = meshgrid(xSample,ySample);

c.X = sum(sum(img .* DistanceX));
c.Y = sum(sum(img .* DistanceY));

end
