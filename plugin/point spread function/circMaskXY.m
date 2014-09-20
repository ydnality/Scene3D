% Create a 2d binary mask for value < reference (1 by default)
%for cartesian coordinate

function Mask=circMaskXY(x,y,varargin)

%INPUT
%R: radius [NxM]
%varargin  {1}: threshold

%OUTPUT
%Mask: binary mask [NxM]


if nargin>2
    th=varargin{1};
else
    th=1; %useful for normalized radius
end

%% Meshgrid the input
[X,Y]=meshgrid(x,y);

R=sqrt(X.*X+Y.*Y);

%% FIND value less than the threshold
[r,c]=find(R<=th); 

%% OUTPUT
Mask=zeros(size(R,1),size(R,2));

for ni=1:length(r)
    Mask(r(ni),c(ni))=1; %Accepted
end