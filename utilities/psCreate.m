function points = psCreate(pX,pY,pZ)
% Create cell array of point sources (field height x depth)
%
%    points = psCreate(pX,pY,pZ)
%
% The point positions are specified.  Nothing is specified about their
% wavelengths.  This is specified elsewhere.
%
% The image side has negative numbers, so pZ < 0
%
% Examples:
%   pZ = -100:20:-40
%   points = psCreate([],[],pZ);
%
%   pX = 0:10:30
%   points = psCreate(pX,[],pZ);
%
% See also: p_re
%
% BW Copyright Vistasoft Team, 2014

if ieNotDefined('pX'), pX = 0; end
if ieNotDefined('pY'), pY = 0; end
if ieNotDefined('pZ'), pZ = -50; end

nFH    = length(pX) * length(pY);
nDepth = length(pZ);

% At each depth, we have a set of 
[X,Y] = meshgrid(pX,pY);
X = X(:); Y = Y(:);
% I don't understand this (BW)
% normalizingZ = min(pZ)*10;
normalizingZ = min(pZ);

% Make the cell array
points = cell(nFH,nDepth);
for ii=1:nFH
    for dd = 1:nDepth
        points{ii,dd} = [X(ii)*  pZ(dd)/normalizingZ, Y(ii), pZ(dd)]; 
    end
end

end

