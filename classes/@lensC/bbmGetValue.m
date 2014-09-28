
function val=bbmGetValue(obj,fileType)

% Get the corresponding field value for th given Black Box Model 
%
%  function val=bbmGetValue(BBoxModel,fileType)
%
%
%INPUT

% 
%fileType: specify which field {'all';'focallength';'focalradius';'imagefocalpoint';'imageprincipalpoint';'imagenodalpoint';'objectfocalpoint';'objectprincipalpoint';'objectnodalpoint';'abcd'}
%
% OUTPUT
%   val : selected value
%
% MP Vistasoft 2014

% Get Black Box equivalente Model
BBoxModel=obj.BBoxModel;
%BBoxModel: struct
%              .focal.length: focal length; 
%              .focal.radius: focal plane radius;
%              .imSpace.focalPoint: focal point in image space
%              .imSpace.principalPoint: principal point in image space
%              .imSpace.nodalPoint: principal point in image space
%              .obSpace.focalPoint: focal point in object space
%              .obSpace.principalPoint: principal point in object space
%              .obSpace.nodalPoint: principal point in object space
%              .abcdMatrix: abcd equivalent matrix; 

switch fileType
    case {'focal'}
        val.length=BBoxModel.focal.length; %focal length
        val.radius=BBoxModel.focal.radius; %focal plane radius
        val.imPoint=BBoxModel.imSpace.focalPoint; %focal point in image space
        val.obPoint=BBoxModel.obSpace.focalPoint; %focal point in object space
    case {'effectivefocallength';'efl'}
        val=BBoxModel.focal.length; %focal length
    case {'focalradius'}
        val=BBoxModel.focal.radius; %focal plane radius
    case {'focalpoint'}
        val.imPoint=BBoxModel.imSpace.focalPoint; %focal point in image space
        val.obPoint=BBoxModel.obSpace.focalPoint; %focal point in object space
    case {'imagefocalpoint'}
        val=BBoxModel.imSpace.focalPoint; %focal point in image space
    case {'objectfocalpoint'}
        val=BBoxModel.imSpace.focalPoint; %focal point in object space
    case {'principalpoint'}
        val.imSpace=BBoxModel.imSpace.principalPoint; %principal point in image space
        val.obSpace=BBoxModel.obSpace.principalPoint; %principal point in object space
    case {'imageprincipalpoint'}
        val=BBoxModel.imSpace.principalPoint; %principal point in image space
    case {'objectprincipalpoint';'objprincipalpoint'}
        val=BBoxModel.obSpace.principalPoint; %principal point in object space
    case {'nodalpoint'}
        val.imSpace=BBoxModel.imSpace.nodalPoint; %principal point in image space
        val.obSpace=BBoxModel.obSpace.nodalPoint; %principal point in image space
    case {'imagenodalpoint'}
        val=BBoxModel.imSpace.nodalPoint; %principal point in image space
    case {'objectnodalpoint';'objnodalpoint'}
        val=BBoxModel.imSpace.nodalPoint; %principal point in image space
    case {'abcd';'abcdmatrix';'abcdMatrix'}
        val=BBoxModel.abcdMatrix; %abcd Matrix 
    case {'all'}
        val=BBoxModel; %all struct
    otherwise
        error (['Not valid: ',fileType,' as field of Black Box Model'])
end
