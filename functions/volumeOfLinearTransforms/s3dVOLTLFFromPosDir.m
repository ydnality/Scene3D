function [ LF] = s3dVOLTLFFromPosDir(XYPos, direction )
% Returns a light-field(LF) matrix from a given position and direction.  
%
% [ LF] = s3dVOLTLFFromPosDir(XYPos, direction )
%
% Inputs
% XYPos: a 2xn matrix containing the position of the rays intersecting a
% known plane. 
% direction: a 2xn matrix containing the direction of the rays.  These are
% the first 2 directions of a unit vector.  
%
% Outputs
% LF: A LF matrix is in the following form: 4xn, where n is the number of
% rays. The first row signifies the X position, the second row the Y
% position, the 3rd row the X direction, the 4th row the Y direction.

 LF = [XYPos(1,:);
        XYPos(2,:);
        direction(1, :);
        direction(2,:)];
end

