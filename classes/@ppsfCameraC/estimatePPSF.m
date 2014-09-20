function ppsfReturn = estimatePPSF(obj,nLines, jitterFlag)
% Calculate the origin and direction of the exiting rays
%
%    ppsfReturn = estimatePPSF(obj,nLines, jitterFlag)
%
% nLines is the number of lines to draw on the diagram.
% For no diagram set nLines to 0 (false).  This is the default.
%
% Example:
%    ppsfCamera.estimatePPSF(nLines)
%
% AL, Vistasoft Copyright 2014

if ieNotDefined('nLines'), nLines = false; end
if ieNotDefined('jitterFlag'), jitterFlag = false; end

disp('-----trace source to lens-----');
tic
ppsfObjectFlag = true;
obj.ppsfRays = obj.lens.rtSourceToEntrance(obj.pointSource, ppsfObjectFlag, jitterFlag);
toc

%duplicate the existing rays, and creates one for each
%wavelength
disp('-----expand wavelenghts-----');
tic
obj.ppsfRays.expandWavelengths(obj.lens.wave);
toc

%lens intersection and raytrace
disp('-----rays trace through lens-----');
tic
obj.lens.rtThroughLens(obj.ppsfRays,nLines);
toc

%project rays onto the z = 0 plane for a proper light field
obj.ppsfRays.projectOnPlane(0);
obj.ppsfRays.pointSourceLocation = obj.pointSource;
ppsfReturn = obj.ppsfRays;

end