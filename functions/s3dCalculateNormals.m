function [normalMap scaledNormalMap z] = s3dCalculateNormals(inputDepthMap, fieldOfView)
% Estimate the surface normals using the depth map
% [normalMap scaledNormalMap] = s3dCalculateNormals(inputDepthMap, fieldOfView)
%
% fieldOfView: horizontal field of view of camera used for depth map(in degrees)
% We perform surface normal estimation by taking the cross product of 2
% perpendicular vectors on the surface.  This vector is averaged amongst
% the 4 pairs of vectors around 1 point of the surface.

% test case - use GT depth map to calculate depth map
% this allows us to decouple the effect of a bad depth map calculation
% method with the 2flashdepth algorithm

% calculate normal vectors using averaged cross products
recordedLateralDistance = zeros(size(inputDepthMap));
normalMap = zeros([size(inputDepthMap,1) size(inputDepthMap,2) 3]);

%convert input depth map to z, instead of depth (there IS a difference!)

%oneD is the assumed D when all points are assumed to be on a plane at z =
%1.  We use similar triangles to figure out z.  1/fakeD = z/d. Here we are
%assuming radial symmetry

z1MaxX = tan(fieldOfView/2 * pi/180);
%assumes that the field of view is for the x dir
z1MaxY = z1MaxX * size(inputDepthMap,1)/size(inputDepthMap,2);  
z1XValues = linspace(-z1MaxX, z1MaxX, size(inputDepthMap,2));
z1YValues = linspace(-z1MaxX, z1MaxY, size(inputDepthMap,1));
[z1XGrid z1YGrid] = meshgrid(z1XValues, z1YValues);
diagDistance = sqrt(z1XGrid.^2 + z1YGrid.^2);  %distance from point on plane to middle of z=1 plane
z1D = sqrt(diagDistance.^2 + 1);

z = inputDepthMap./z1D;
figure; imagesc(z); title('z values');


for i = 2:(size(z, 2) - 1)
    for j = 2:(size(z,1) -1)
        aRelief = z(j - 1,i) - z(j,i);
        bRelief = z(j, i + 1) - z(j,i);
        cRelief = z(j + 1,i) - z(j,i);
        dRelief = z(j,i - 1) - z(j,i);
        
        %lateralDistance = inputDepthMap(j,i) * pixelUnit/sensorDistance;
        %%this does not make a lot of sense
        
%         lateralDistance = inputDepthMap(j,i) * tan(fieldOfView/size(ratioImage,2) * pi/180);
        lateralDistance = z(j,i) * tan(fieldOfView/size(inputDepthMap,2) * pi/180);
        recordedLateralDistance(j,i) = lateralDistance;
        
        aVector = [0 lateralDistance -aRelief];
        bVector = [lateralDistance 0 -bRelief];
        cVector = [0 -lateralDistance -cRelief];
        dVector = [-lateralDistance 0 -dRelief];
        
        adNormal = cross(aVector, dVector);
        adNormal = adNormal./norm(adNormal);
        
        dcNormal = cross(dVector, cVector);
        dcNormal = dcNormal./norm(dcNormal);
        
        cbNormal = cross(cVector, bVector);
        cbNormal = cbNormal./norm(cbNormal);
        
        baNormal = cross(bVector, aVector);
        baNormal = baNormal./norm(baNormal);
        
        averageNormal = (adNormal + dcNormal + cbNormal + baNormal);
        averageNormal = averageNormal./norm(averageNormal);
        
        normalMap(j, i,:) = reshape(averageNormal, [1 1 3]);
    end
end
%normalize the normal vectors;
normalMap = normvec(normalMap, 'p', 2, 'dim', 3);
scaledNormalMap = normalMap./2 + .5;

