function [csBasis,cSignals,sValues] = lmColorSignalBasis(lights,surfaces,nBases)
%
%  [csBasis,cSignals,sValues] = lmColorSignalBasis(lights,surfaces,nBases)
%
% Calculate the color signal basis from light SPDs and surface reflectances
% (in the columns). If requested, return all of the color signals, as well.
%
% Example:
%  lights = lmLookupSignals({'FL7','Vivitar'}, 'illuminants',[400:10:700],1);
%  surfaces = lmLookupSignals(surfaceList, 'surfaces', wavelength);
%  [csBasis,cs,sValues] = lmColorSignalBasis(lights,surfaces,nBases)
%
%

nLights = size(lights,2);
nWave   = size(lights,1);
if size(surfaces,1) ~= nWave
    error('Surface and lights must have the same wavelength representation.');
end

if find(isnan(lights)),   error('NaN in the lights'); end
if find(isnan(surfaces)), error('NaN in the surfaces'); end

%  Multiply the surfaces by each of the lights and combine them into a
%  single color signals matrix.
cSignals = [];
for ii=1:nLights
    cSignals = [cSignals, diag(lights(:,ii))*surfaces];
end

if ~exist('nBases','var'),  nBases = []; end
[csBasis,sValues] = lmComputeBases(cSignals,0,nBases);

% plot(csBasis); grid on

return;
