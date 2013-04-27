function oi = pbrt2oi(fname)
% %Convert pbrt multispectral irradiance data into an ISET scene
% %
% %   oi = pbrt2oi(fname)
% %
% % Example
% %  fname = 'teapot-subsurface-small.dat';
% %  oi = pbrt2oi(fname);
% %  vcAddAndSelectObject(oi); oiWindow;
% %
% % (c) Stanford VISTA Team 2012
% 
% % TODO
% %   Add in lens and other information to the pbrt file and incorporate it.
% %
% 
% if ieNotDefined('fname'), error('File name required.'); end
% 
% %% Read header
% fID = fopen(fname,'r','l');
% [A, count] = fscanf(fID,'%d %d %d\n',[3 Inf]);
% if count ~= 3
%     fclose(fID);
%     error('Misread header.'); 
% end
% 
% %% Read photons
% [photons, count] = fread(fID,prod(A),'double');
% if count ~= prod(A)
%     fclose(fID);
%     error('Misread photons');
% end
% 
% % Close file
% fclose(fID);
% 
% 
% %% Check dynamic range of photon count.
% mn = min(photons); mx = max(photons);
% if mx/mn > 2^16
%     warndlg('Image dynamic range is very high');
% end
% 
% %% Set the OI data in the OI structure
% oi = oiCreate('uniform d65');
% % Set the optics parameters from somewhere
% % optics = oiGet(oi,'optics');
% photons = reshape(photons,A(1),A(2),A(3));
% 
% oi = oiSet(oi,'cphotons',photons);
% 
% oi = oiSet(oi,'name',fname);
% oi = oiSet(oi,'h fov',10);
% oi = oiSet(oi,'depth map',[]);
% 
% illum = oiCalculateIlluminance(oi);
% oi = oiSet(oi,'illuminance',illum);
% 
% % Set mean illuminance to 10 lux
% oi = oiAdjustIlluminance(oi,10);


%Convert pbrt multispectral irradiance data into an ISET scene
%
%   oi = pbrt2oi(fname)
%
%

if ieNotDefined('fname'), error('File name required.'); end

%open file
fID = fopen(fname,'r','l');

%load size information
[A, cnt] = fscanf(fID,'%d %d %d\n',[3 1]);

%load lens and field of view and information
[FOV, cnt2] = fscanf(fID,'%f %f %f\n',[3 1]);

if (~isempty(FOV))
    focalLength = FOV(1)
    aperture = FOV(2)
    fiedOfView = FOV(3)
else
    'no lens information!!'
end


%Load the stored photons produced by AL's pbrt code
photons = fread(fID,prod(A),'double');

%A(2) = rows, A(1) = columns
photons = reshape(photons,A(2),A(1),A(3));
fclose(fID);

% Set the OI data
oi = oiCreate;
oi = initDefaultSpectrum(oi);
% Set the optics parameters from somewhere
% optics = oiGet(oi,'optics');

oi = oiSet(oi,'cphotons',photons);

%temporarily disable depth map reading
% [path,name,ext] = fileparts(fname); 
% dMapFile = [name '_dm_DM.dat']; 
% oi.depthMap = s3dReadDepthMapFile(dMapFile);
%oi.depthMap = imresize(scene.depthMap, [sceneGet(scene, 'rows') sceneGet(scene, 'cols')]);

return
