function [LFout, LFmid, LFin] = s3dLightField(pointSource, lens)
% Create a light field object for the input (point source) and lens
%
% LF = s3dLightField(pointSource, lens)
%
% pointSOurce:
% lens
% film
%
% LF:  Light field object
%
% Example:
%  lens = lensC; pointSource = [0 1.7 -103];
%  [LFout, LFmid, LFin] = s3dLightField(pointSource, lens);
%
% See Also:
%
% AL, VISTASOFT, 2014

%% ray trace and save ppsf - Not sure camera should have pointSources

% Use the multi element lens and film and a point source.  Combine into
% a camera that calculates the point spread function.
film = filmC;

ppsfCamera = ppsfCameraC('lens', lens, 'film', film, 'pointSource', pointSource);

% Plenoptic point spread calculated with Snell's Law
ppsf = ppsfCamera.estimatePPSF(0, true);  %0 debug lines; jitter set to true


%% Calculate light fields
    
LFin  = ppsf.LF('in');
LFmid = ppsf.LF('middle');
LFout = ppsf.LF('out');

end
