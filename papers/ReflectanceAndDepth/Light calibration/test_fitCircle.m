close all;
clear variabels;
clc;

xc = -30;
yc = 10;
r = [47 36 18];

xCoord = [];
yCoord = [];
label = [];

for i=1:length(r)
    angles = 1:90;
    nAngles = length(angles);

    xCoord = [xCoord r(i)*cosd(angles)  + xc + randn(1,nAngles)];
    yCoord = [yCoord r(i)*sind(angles) + yc + randn(1,nAngles)];
    label = [label i*ones(nAngles,1)];
end

figure; 
hold on; grid on; box on;
plot(xCoord,yCoord,'x');
plot(xc,yc,'gs');
axis square; axis equal;


[xEst, yEst, rEst] = fitConcentricCircles(xCoord',yCoord',label);


plot(xEst,yEst,'ro')
for i=1:length(rEst)
    rectangle('Position',[xEst-rEst(i) yEst-rEst(i) 2*rEst(i) 2*rEst(i)],'Curvature',[1 1]);
end