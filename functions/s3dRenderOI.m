function oi = s3dRenderOI(inputPbrt, oiName, dockerFlag, noScale)
%
%    oi= s3dRenderOI(inputPbrt, oiName, dockerFlag, noScale)
%
% inputPbrt:   pbrt object or a full path to a file that describes the pbrt
%              data.
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

if (ieNotDefined('oiName')),       oiName = 'default'; end
if (ieNotDefined('dockerFlag')),   dockerFlag = true;  end
if ieNotDefined('noScale'),        noScale = false; end

oiFlag = true;

% Render using the docker calculation
if (isa(inputPbrt, 'pbrtObject'))
    oiPbrt = pbrtObject;
    oiPbrt.makeDeepCopy(inputPbrt);
else
    oiPbrt = inputPbrt;
end

oi = s3dRenderScene(oiPbrt, oiName, noScale, dockerFlag, oiFlag);
% vcAddObject(oi); oiWindow;
 
end