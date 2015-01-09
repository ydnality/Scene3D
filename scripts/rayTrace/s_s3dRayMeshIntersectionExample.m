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


%% Now try an example that is more realistic for our applications (simple depth map with one ray)

% Simple depth map
[x, y] = meshgrid(1:100,1:100);
zMap = ones(100) * -50;
zMap(25:75, 25:75) = -25;
vcNewGraphWin; mesh(zMap);

faces = delaunay(x,y);       % net list for triangles
orig  = [20 20 -50];         % ray's origin
dir   = [50 50 0] - orig;    % ray's direction

% calculate the intersection from a ray to the depth map

vertices = [x(:) y(:) zMap(:)];
vert1 = vertices(faces(:,1), :);  %get vertex coordinates of triangles
vert2 = vertices(faces(:,2), :);
vert3 = vertices(faces(:,3), :);

[intersect,~,~,~,xcoor] = TriangleRayIntersection(orig, dir, ...
  vert1, vert2, vert3, 'lineType' , 'line');

%print intersection information
fprintf('Number of: faces=%i, points=%i, intresections=%i; time=%f sec\n', ...
  size(faces,1), size(vertices,1), sum(intersect), toc);

% Plot the intersection
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



%% Now try a depth map with multiple rays

% Simple depth map
[x, y] = meshgrid(1:100,1:100);
zMap = ones(100) * -50;
zMap(25:75, 25:75) = -25;
vcNewGraphWin; mesh(zMap);

faces = delaunay(x,y);       % net list for triangles
orig  = [20 20 -50];         % ray's origin
% [xOrig yOrig zOrig] = meshgrid(15:1:25, 15:1:25, -50);
% xOrig = xOrig(:);
% yOrig = yOrig(:);
% zOrig = zOrig(:);
% orig = [xOrig yOrig zOrig];
%dir   = repmat([50 50 0], [size(xOrig,1), 1]) - orig;    % ray's direction

%sample a square aperture centered at the center of the lens
[xDir yDir zDir] = meshgrid(45:55, 45:55, 0);
xDir = xDir(:);
yDir = yDir(:);
zDir = zDir(:);

dir = [xDir yDir zDir];
orig = repmat(orig, [size(xDir, 1), 1]);
dir = dir - orig;

% calculate the intersection from a ray to the depth map

vertices = [x(:) y(:) zMap(:)];
vert1 = vertices(faces(:,1), :);  %get vertex coordinates of triangles
vert2 = vertices(faces(:,2), :);
vert3 = vertices(faces(:,3), :);

%clone so that orig, dir, and vertn are Nx3 in dimensions
% orig = [orig; repmat(orig(size(orig,1), :), [size(vert1,1) - size(orig,1) 1])];
% dir = [dir; repmat(dir(size(dir,1), :), [size(vert1,1) - size(dir,1) 1])];
% Andy: this doesn't really work

% what we really need is to clone both vertices, and orig/dir, so that each
% orig/dir entry sees ALL triangles as potential intersection locations,
% not just one.  We do this by repeating the rays, and for each ray have a
% potential intersection with each one of the triangles.


%total potential intersections = size(rays) * size(vertices)
%so rays are duplicated #vertices times, and vertices are duplicated #ray
%times

origExp = repmat(orig, [size(vert1,1) 1]);  %repeat origin so can potentially intersect with all triangles
dirExp = repmat(dir, [size(vert1,1) 1]);    %cycles once before repeating

vert1Exp = repmat(vert1(:), [1 size(orig, 1)])';    %repeats each triangle size(orig) times
vert1Exp = reshape(vert1Exp, size(dirExp,1), [] );     %ex.  1 1 2 2 3 3 4 4 

vert2Exp = repmat(vert2(:), [1 size(orig, 1)])'; 
vert2Exp = reshape(vert2Exp, size(dirExp,1), [] );  

vert3Exp = repmat(vert3(:), [1 size(orig, 1)])'; 
vert3Exp = reshape(vert3Exp, size(dirExp,1), [] );  


[intersect,~,~,~,xcoor] = TriangleRayIntersection(origExp, dirExp, ...
  vert1Exp, vert2Exp, vert3Exp);%, 'lineType' , 'line');

%print intersection information
fprintf('Number of: faces=%i, points=%i, intersections=%i; time=%f sec\n', ...
  size(faces,1), size(vertices,1), sum(intersect), toc);


% Plot the intersection
%plot the surface, with intersecting triangles in red
vcNewGraphWin;
%trisurf(faces, x,y,zMap, intersect*1.0, 'FaceAlpha', .5);
trisurf(faces, x,y,zMap,  'FaceAlpha', .5);
hold on;
%plot intersection coordinates
scatter3(xcoor(intersect,1), xcoor(intersect,2), xcoor(intersect,3), 100, 'b', 'o', 'filled')
%plot ray

% for i = 1:size(orig,1)
% hold on;
% line('XData',orig(i, 1)+[0 dir(i,1)],'YData',orig(i,2)+[0 dir(i,2)],'ZData',...
%   orig(3)+[0 dir(3)],'Color','r','LineWidth',3)
% end




%% Now try ray for a hemisphere, kind of like a lens (brian's example).

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

