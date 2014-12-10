%% Example of ray(s) - mesh intersection 
%
% A lot of code is borrowed from ray/triangle intersection toolbox
% developed by Jarek Tuszynski* (jaroslaw.w.tuszynski@leidos.com)
%
% This script illustrates how to make a ray, a mesh, and run the ray/mesh
% intersection code in the external toolbox. 
%
% The toolbox is added in external/rayMesh...
%
% AL Vistasoft Team 2014

%% Create an example mesh and intersects a ray

n=20;
[x,y] = meshgrid(1:n,1:n);    % create 2D mesh of points
faces = delaunay(x,y);        % triangulate the mesh using Delaunay 
z     = peaks(n);             % sample function defined on a grid of the same dimenision
vertices = [x(:) y(:) z(:)];  % vertices stored as Nx3 matrix


% For this toolbox the direction is the end point minus the origin.
% So, the endpoint is dir + origin
orig  = [0.25*n 0 2];         % ray's origin
dir   = [0.5 *n n 0];         % ray's direction
vert1 = vertices(faces(:,1),:);
vert2 = vertices(faces(:,2),:);
vert3 = vertices(faces(:,3),:);

% Fast
tic;
intersect = TriangleRayIntersection(orig, dir, vert1, vert2, vert3);
fprintf('Number of: faces=%i, points=%i, intresections=%i; time=%f sec\n', ...
  size(faces,1), size(vertices,1), sum(intersect), toc);

%% Display the results
% Surface in blue, line in light read and intersected
% triangles in dark red*
figure(1); clf;
trisurf(faces,x,y,z, intersect*1.0,'FaceAlpha', 0.9)
hold on;
line('XData',orig(1)+[0 dir(1)],'YData',orig(2)+[0 dir(2)],'ZData',...
  orig(3)+[0 dir(3)],'Color','r','LineWidth',3)
set(gca, 'CameraPosition', [106.2478  -35.9079  136.4875])
%set(gco,'EdgeColor','none');


%% Example with many rays and many triangles (many faces / many rays type problem)
% So far all examples were of a single ray and many triangles. However
% one can as well have one triangle and many rays, or many rays and many 
% triangles. 
%
% Example below calculates intersections between faces and rays 
% going through the center of each face. Since each intersection is in the 
% same relative point t, u and v returned are very similar. Plot shows 
% intersection points

faceCenter = (vert1+vert2+vert3)/3;
Orig  = repmat(orig,size(vert1,1),1); % Clone it until the same size as vert1
[intersect, t, u, v, xcoor] = TriangleRayIntersection(Orig, ...
  2*(faceCenter-Orig), vert1, vert2, vert3);
fprintf('Number of: faces=%i, intresections=%i\n', size(faces,1), sum(intersect));
fprintf('mean t=%f+-%f\n', mean(t), std(t));
fprintf('mean u=%f+-%f\n', mean(u), std(u));
fprintf('mean v=%f+-%f\n', mean(v), std(v));

figure(1); clf;
plot3(xcoor(:,1), xcoor(:,2), xcoor(:,3), 'o')


%% Now try an example that is more realistic for our applications

% Simple depth map
[x, y] = meshgrid(1:100,1:100);
zMap = ones(100) * -50;
zMap(25:75, 25:75) = -25;
vcNewGraphWin; mesh(zMap);

faces = delaunay(x,y);       % net list for triangles
orig  = [20 20 -50];         % ray's origin
dir   = [50 50 0] - orig;    % ray's direction

%% calculate the intersection from a ray to the depth map

vertices = [x(:) y(:) zMap(:)];
vert1 = vertices(faces(:,1), :);  %get vertex coordinates of triangles
vert2 = vertices(faces(:,2), :);
vert3 = vertices(faces(:,3), :);

[intersect,~,~,~,xcoor] = TriangleRayIntersection(orig, dir, ...
  vert1, vert2, vert3, 'lineType' , 'line');

%print intersection information
fprintf('Number of: faces=%i, points=%i, intresections=%i; time=%f sec\n', ...
  size(faces,1), size(vertices,1), sum(intersect), toc);

%% Plot the intersection
%plot the surface, with intersecting triangles in red
vcNewGraphWin;
trisurf(faces, x,y,zMap, intersect*1.0, 'FaceAlpha', .5);
hold on;
%plot intersection coordinates
scatter3(xcoor(intersect,1), xcoor(intersect,2), xcoor(intersect,3), 100, 'b', 'o', 'filled')
%plot ray
hold on;
line('XData',orig(1)+[0 dir(1)],'YData',orig(2)+[0 dir(2)],'ZData',...
  orig(3)+[0 dir(3)],'Color','r','LineWidth',3)

%% Now try multiple rays for a hemisphere, kind of like a lens.

% Faces and vertices of the lens
rCurve = 1.8; Ap = 1; nSamples = 100;
[f,v] = mshSphericalLens(rCurve,Ap,nSamples);

vcNewGraphWin;
trisurf(f,v(:,1),v(:,2),v(:,3), 'FaceAlpha', .5); axis equal

% The ray
orig = [0 0 -3];
dir  = [.5 .5 1] - orig;

% vN: the Nth entry of the triangle vertices
% So vN is nFaces x 3
%
v1 = v(f(:,1), :);  %get vertex coordinates of triangles
v2 = v(f(:,2), :);
v3 = v(f(:,3), :);
[intersect,~,~,~,xcoor] = TriangleRayIntersection(orig, dir, ...
  v1,v2,v3, 'lineType' , 'line');

vcNewGraphWin;
trisurf(f,v(:,1),v(:,2),v(:,3), intersect*1.0, 'FaceAlpha', .5);
hold on; axis equal

%plot intersection coordinates
scatter3(xcoor(intersect,1), xcoor(intersect,2), xcoor(intersect,3), 100, 'b', 'o', 'filled')

%plot ray
hold on;
line('XData',orig(1)+[0 dir(1)],'YData',orig(2)+[0 dir(2)],'ZData',...
  orig(3)+[0 dir(3)],'Color','r','LineWidth',3)

axis equal



%% END

