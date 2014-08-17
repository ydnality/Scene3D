% Create a 2d binary mask for value < reference (1 by default)

function Mask=circMask(R,varargin)

%INPUT
%R: radius [NxM]
%varargin  {1}: threshold

%OUTPUT
%Mask: binary mask [NxM]


if nargin>1
    th=varargin{1};
else
    th=1; %useful for normalized radius
end


%% FIND value less than the threshold
[r,c]=find(R<=th); 

%% OUTPUT
Mask=zeros(size(R,1),size(R,2));

for ni=1:length(r)
    Mask(r(ni),c(ni))=1; %Accepted
end
