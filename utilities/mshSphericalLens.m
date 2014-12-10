function [faces,vertices] = mshSphericalLens(RoC,Ap,n)
% Make a hemisphere with a specific radius of curvature and aperture
%
%
% RoC - Radius of curvature
% Ap  - Aperture radius
% n   - Number of sample points on the sphere
%
% Triangular faces and vertices
%  faces
%  v1,v2,v3:  
%
% Example:
%  [f,v] = mshSphericalLens(1.5,0.5,100);
%  vcNewGraphWin; trisurf(f, v(:,1),v(:,2),v(:,3), 'FaceAlpha', .5); axis equal; hold on;
%
% AL Vistasoft Copyright 2015

[X,Y,Z] = sphere(n);

% Find half the sphere
lst = (Z + eps < 0);
X = X(lst); Y = Y(lst); Z = Z(lst);

% Scale for the radius of curvature
X = X(:)*RoC; Y = Y(:)*RoC; Z = Z(:)*RoC;

% Find the vertices inside the aperture
pos = sqrt(X(:).^2 + Y(:).^2);
lst = (pos <= Ap);
X = X(lst); Y = Y(lst); Z = Z(lst);

% Set up the vertices
vertices = [X(:),Y(:),Z(:)];

% Make the triangular faces from the vertices
faces = delaunay(vertices(:,[1 2]));

end

