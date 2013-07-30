function signals = lmLookupSignals(signalList, dbName, wavelength, normalize)
%
%  signals = lmLookupSignals(signalList, dbName, wavelength)
%
%Author: FX, BW
%Purpose:
%    Caclulate the spectral base functions for color signals based on the
%    assumption of surfaces and lighting statistics
%
% Example:
%   lights = lmLookupSignals({'FL7','Vivitar'}, 'illuminants',[400:10:700],1);
%   lights = lmLookupSignals({'FL7','Vivitar'}, 'illuminants',[400:10:700],1);

if ~exist('wavelength','var'),  wavelength = 400:10:700; end
if ~exist('normalize','var'),  normalize = 0; end

% Read in all the signals from the data base
signals = [];
for ii=1:length(signalList)
    signals  =  [signals,vcReadSpectra(signalList{ii},wavelength)];
end

if normalize
    % normalize the total energy for each light
    signals = signals/diag(sum(signals));
end

return;


