



function  [lens]  = lensBBoxModelOLD(lens,OptSyst)


% Get the parameters from the optical system structure to build an
% equivalent Black Box Model of the lens.



% INPUT
% lens:
% OptSyst:



%OUPUT
% lens:   added with an equivalent Black Box Model



%% SET OUTPUT - this could be the optConvert routine
        
% Get 'new' origin for optical axis 
z0 = OptSyst.cardPoints.lastVertex;

lens.BBoxModel.focal.length = OptSyst.cardPoints.fi; %focal lenght of the system
lens.BBoxModel.focal.Radius = OptSyst.Petzval.radius; % radius of curvature of focal plane
% The result is an equivalent lens formulation.
% The new formulation is specified using principal points and nodal
% points and focal points
%
% pLensC = lensPrincipal(lens);

lens.BBoxModel.ImSpace.focalPoint=OptSyst.cardPoints.dFi;     %Focal point in the image space
lens.BBoxModel.ImSpace.principalPoint=OptSyst.cardPoints.dHi; % Principal point in the image space
lens.BBoxModel.ImSpace.nodalPoint=OptSyst.cardPoints.dNi;     % Nodal point in the image space
lens.BBoxModel.ObjSpace.focalPoint=OptSyst.cardPoints.dFo-z0; %Focal point in the object space
lens.BBoxModel.ObjSpace.principalPoint=OptSyst.cardPoints.dHo-z0; % Principal point in the object space
lens.BBoxModel.ObjSpace.nodalPoint=OptSyst.cardPoints.dNo-z0; % Nodal point in the object space

lens.BBoxModel.abcdMatrix = OptSyst.matrix.abcd; % The 4 coefficients of the ABCD matrix of the overall system


