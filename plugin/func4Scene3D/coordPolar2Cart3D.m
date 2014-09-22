

function  [x,y,zOUT]=coordPolar2Cart3D(ro,theta,z)

%Convert point coordinate from  polar to Cartesian (x,y,z)  in object/image plane(ro,theta,z)

% INPUT:
% ro= height eccentricity
% theta= angular eccentricity
% z= distance along optical axis 

% OUTPUT:
% x
% y
% z






%% X coord
x=ro.*cos(theta);
%% Y coord
y=ro.*sin(theta);
%% Z coord
zOUT=z;