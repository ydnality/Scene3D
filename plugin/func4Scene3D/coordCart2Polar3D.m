

function  [ro,theta,zOUT]=coordCart2Polar3D(x,y,z)

%Convert point coordinate from Cartesian (x,y,z) to polar in object/image plane(ro,theta,z)

% INPUT:
% x
% y
% z

% OUTPUT:
% ro= height eccentricity
% theta= angular eccentricity
% z= distance along optical axis 


%% HEIGHT ECCENTRICITY
ro =sqrt(x.^2+y.^2); %distance of the point source to optical axis

%% ANGULAR ECCENTRICITY
if not(x==0) && not(y==0)
    theta=atan(y/x); %angle subtended by the point source and the x-axis in the object plane
else
    if (y==0)
        theta=0;
    else
        ps_angle=pi/2;
    end
end
        
%% DISTANCE along THE OPTICAL AXIS
zOUT=z;