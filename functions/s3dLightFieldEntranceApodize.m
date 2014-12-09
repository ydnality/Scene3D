function [LFin] = s3dLightFieldEntranceApodize(pointSource, lens, depthMap)
% Create a light field object for the input (point source) and lens
%
%
% General idea of this function:
%
% This function performs the same function as s3dLightFieldEntrance, except
% that it allows for the proper rendering of partially occluded regions,
% assuming that we have the proper information in those regions.

% For each ray that is shot at the aperture, we will project that ray onto
% the x-y plane.  For sampled places on this ray, we will compare and see
% if the ray is greater than or less than the depth of the scene.  If it is
% further from the scene than the depth, then that particular ray is
% occluded.
%
% The amount of sampling is important because we can potentially miss
% crucial parts of the scene if we sample too coursely.
%
% LF = s3dLightFieldApodize(pointSource, lens, depth)
%
% pointSOurce:
% lens
% film
%
% LF:  Light field object
%
% Example:
%  lens = lensC; pointSource = [0 1.7 -103];
%  [LFout, LFmid, LFin] = s3dLightField(pointSource, lens);
%
% See Also:
%
% AL, VISTASOFT, 2014

%% ray trace and save ppsf - Not sure camera should have pointSources

% Use the multi element lens and film and a point source.  Combine into
% a camera that calculates the point spread function.
film = pbrtFilmC;

ppsfCamera = ppsfCameraC('lens', lens, 'film', film, 'pointSource', pointSource);



% Plenoptic point spread calculated with Snell's Law
% change this function here. 

% we need a ppsfCamera.traceToEntrance(0, true, z);  %this is depth aware

% before we obtain z, we need to convert the depth map, from depth map
% format to z format. How do we properly interpolate values for z where
% there is no value???

%embrace it and keep in point cloud format?


%load seamount
%tri = delaunay(x,y);
%trisurf(tri,x,y,z);

%TriScatteredInterp


%Create a data set:
x = rand(100,1)*4-2;
y = rand(100,1)*4-2;
z = x.*exp(-x.^2-y.^2);

%Construct the interpolant:
F = TriScatteredInterp(x,y,z);

%Evaluate the interpolant at the locations (qx, qy). The corresponding value at these locations is qz:
ti = -2:.25:2;
[qx,qy] = meshgrid(ti,ti);
qz = F(qx,qy);
mesh(qx,qy,qz);
hold on;
plot3(x,y,z,'o');

ppsf = ppsfCamera.traceToEntrance(0, true);  %0 debug lines; jitter set to true


%% Calculate light fields
    
LFin  = ppsf.LF('in');

end
