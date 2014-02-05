function pbrtPath = pbrtRootPath
% Return the path to the root directory of pbrt on this system.
%
%   pbrtPath = pbrtRootPath
%
% This should be a which pbrt call, really.
% I am also not sure why it is here.  Why don't we just ask users to put
% pbrt on their path sso unix('pbrt') runs?
%
%s

if exist('/home/andy/','dir')
    pbrtPath = '/home/andy/Dropbox/Scene3D/pbrt-v2-spectral-diffraction/';
elseif exist('/home/ajwandell','dir')
    pbrtPath = '/home/ajwandell/Dropbox/Scene3D/pbrt-v2-spectral-diffraction/';
elseif exist('/celadon','dir')
    pbrtPath = '/celadon/software/pbrt-v2-spectral-diffraction/';
end

return

