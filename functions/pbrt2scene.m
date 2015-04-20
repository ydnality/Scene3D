function scene = pbrt2scene(fname)
%Convert pbrt multispectral irradiance data into an ISET scene
%
%   scene = pbrt2scene(fname)
%
% See also:  pbrt2oi, s3dRenderScene, s3dRenderSceneAndDepthMap
%
% (c) Stanford VISTA Team 2012

%%
if ieNotDefined('fname'), error('File name required.'); end

% Open file
fID = fopen(fname,'r','l');

%load size information
[A, cnt] = fscanf(fID,'%d %d %d\n',[3 1]);

%load lens and field of view and information
[FOV, cnt2] = fscanf(fID,'%f %f %f\n',[3 1]);

% if (~isempty(FOV))
%     focalLength = FOV(1);   %do something with this information in the future
%     aperture = FOV(2);
%     fiedOfView = FOV(3);
% else
%     'no lens information!!'
% end

%%Load the stored photons produced by the pbrt code
photons = fread(fID,prod(A),'double');

%A(2) = rows, A(1) = columns
photons = reshape(photons,A(2),A(1),A(3));
fclose(fID);

%% Set the scene photon data
scene = sceneCreate;
scene = initDefaultSpectrum(scene);
scene = sceneSet(scene,'photons',photons(:,:,1:31));    %ignore the 32nd wavelength for now

% How do we figure out the field of view and such?  Probablky from  FOV



%% Temporarily disable depth map reading
% [path,name,ext] = fileparts(fname); 
% dMapFile = [name '_dm_DM.dat']; 
% oi.depthMap = s3dReadDepthMapFile(dMapFile);
% oi.depthMap = imresize(scene.depthMap, [sceneGet(scene, 'rows') sceneGet(scene, 'cols')]);

end
