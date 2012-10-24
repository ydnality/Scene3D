function depthCenters = oiCalculateDepthCenters(depthEdges)
% calulates depthCenters, from depth edges (for easier computation later)
% 
% 
% Inputs:
%   depthEdges: Vector of the boundaries of the depth bins. 
% 
% Return values:
%   depthCenters: depth centers (center of bins)
% 
% Examples:  depthCenters = oiCalculateDepthCenters([0 1 2 3 4 5])
% 

if ieNotDefined('depthEdges')
    depthEdges = [min(dMap(:)),max(dMap(:))];
elseif length(depthEdges) == 1
    depthEdges = [min(dMap(:)),depthEdges, max(dMap(:))];
end

depthCenters = zeros(length(depthEdges)-1,1);
for ii=1:(length(depthEdges)-1)
    depthCenters(ii) = depthEdges(ii) + (depthEdges(ii+1) - depthEdges(ii))/2;
end