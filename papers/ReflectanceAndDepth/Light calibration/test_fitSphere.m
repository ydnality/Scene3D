close all;
clear variabels;
clc;

xc = -30;
yc = 10;
zc = 66;
r = 13;
wg = [1 3 5];

xCoord = [];
yCoord = [];
zCoord = [];
wght = [];

for i=1:length(wg)
    [phi, theta] = meshgrid(-10:2:10,0:10:360);
    phi = phi(:)';
    theta = theta(:)';
    nAngles = length(phi);

    xCoord = [xCoord r*wg(i)*sind(theta).*cosd(phi)  + xc + randn(1,nAngles)];
    yCoord = [yCoord r*wg(i)*sind(theta).*sind(phi) + yc + randn(1,nAngles)];
    zCoord = [zCoord r*wg(i)*cosd(theta) + zc + randn(1,nAngles)];
    wght = [wght wg(i)*ones(1,nAngles)];
end

figure; 
hold on; grid on; box on;
scatter3(xCoord,yCoord,zCoord,'x');
axis square; axis equal;


[xEst, yEst, zEst, rEst] = fitConcentricSpheres(xCoord',yCoord',zCoord',wght');
scatter3(xEst, yEst, zEst, 'sm');



plot(xEst,yEst,'ro')
for i=1:length(rEst)
    rectangle('Position',[xEst-rEst(i) yEst-rEst(i) 2*rEst(i) 2*rEst(i)],'Curvature',[1 1]);
end