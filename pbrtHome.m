function [pbrtPath]=pbrtHome
%
%   rootPath = s3dRootPath;
%
% This file must reside in the root directory of the scene3D Project
% folder.  The file is used to identify the directory of PBRT.

error('Obsolete')

%please update this to reflect the proper location of pbrt.
%pbrtPath = '/home/andy/Dropbox/Scene3D/pbrt-v2-spectral-diffraction/';
%pbrtPath = '/home/ydna/pbrt-v2-spectral-diffraction/';

%rootPath=which('pbrtHome');
%[rootPath,fName,ext]=fileparts(rootPath);

if (exist('/home/andy/'))
    pbrtPath = '/home/andy/Dropbox/Scene3D/pbrt-v2-spectral-diffraction/';
else
    pbrtPath = '/home/ajwandell/Dropbox/Scene3D/pbrt-v2-spectral-diffraction/';
end

return;

