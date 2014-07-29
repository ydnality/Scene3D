function pbrtPath = pbrtRootPath
% Return the path to the root directory of pbrt on this system.
%
%   pbrtPath = pbrtRootPath
%
% This should be a which pbrt call, really.
% I am also not sure why it is here.  Why don't we just ask users to put
% pbrt on their path sso unix('pbrt') runs?
%
%
% New solution: Add the pbrt directory to the matlab script.  
% 
% Ex.
%    
% add the following line to the beginning of the matlab script used to
% launch matlab, where <pbrtDirectory> is the directory of the pbrt binary file:
% export PATH=$PATH:<pbrtDirectory>

error('Obsolete')

if exist('/home/andy/','dir')
    pbrtPath = '/home/andy/Dropbox/Scene3D/pbrt-v2-spectral-diffraction/';
elseif exist('/home/ajwandell','dir')
    pbrtPath = '/home/ajwandell/Dropbox/Scene3D/pbrt-v2-spectral-diffraction/';
elseif exist('/celadon','dir')
    pbrtPath = '/celadon/software/pbrt-v2-spectral-diffraction/';
end

return

