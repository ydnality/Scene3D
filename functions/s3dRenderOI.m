function oi = s3dRenderOI(inputPbrt, focalLength, oiName, dockerFlag, noScale)
%
%    oi= s3dRenderOI(inputPbrt, focalLength, oiName, dockerFlag, noScale)
%
% inputPbrt:   pbrt object or a full path to a file that describes the pbrt
%              data.
% focalLength: optics focal length
% 
% sceneName: file name of the pbrt file to render (without the .pbrt
% extension).  sceneName will also be used to name the optical image.
%
% This function renders a PBRT oi using the given pbrt object, then
% returns the data as an optical image. The pbrt file should be previously 
% generated and be in datapath/generatedPbrtFiles/.  The output file is placed
% in datapath/generatedPbrtFiles/tempPBRT/.  The
% proper optics is also placed into this optical image that corresponds to
% the focal length and field of view (FOV).
%
% Example
%  inputPbrt = pbrtObject; dockerFlag = true; noScale = false;
%  oi= s3dRenderOI(inputPbrt, [], 'deleteMe', dockerFlag, noScale)
%  vcAddObject(oi); oiWindow;
% Todo: We are considering copying the pbrt file that was used to generate
% the scene when we save the oi for future use.
% 
% Todo: backwards compatibility to render pbrt files

if (ieNotDefined('oiName')),    oiName = 'default'; end
if (ieNotDefined('dockerFlag')),   dockerFlag = false;  end
if (ieNotDefined('focalLength')),  focalLength = 0.05; end
if ieNotDefined('noScale'),        noScale = false; end

oiFlag = true;
oi = s3dRenderScene(inputPbrt, oiName, noScale, dockerFlag, oiFlag);

% We need the right way to figure out the effective focal length.  If
% inputPbrt is a pbrtObject we can query.  If inputPbrt is a file name, we
% will need to figure it out some other way that we haven't determined yet.
% For now, we hope the user sent it in, or we just make up some B.S. and
% say 50mm.
if ischar(inputPbrt) 
    % We should read the lens file and turn it into a lens object.
    % Then we should run the lens ABCD stuff.
    % Then we should get the equivalent focal length (efl) of the lensC
    % object.

elseif isa(inputPbrt,'pbrtObject')
    % It is  a pbrt Object
    % In some cases we should be able to figure out the effective focal
    % length.
end

oi = oiSet(oi, 'wangular', 5);  %prevent error
oi = s3dFixOi(oi, focalLength);
 
end