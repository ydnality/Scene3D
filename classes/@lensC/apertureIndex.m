function [varargout] = apertureIndex(obj,index)
% Set the surface subtype specify by the index as  'diaphragm' 
% It is required to compute the Black Box Model
%
%
%  In this format, it creates the BBM with n_ob and n_im = 1
%  BBM = lens.set('apertureindex', index)
%    same as 
%  BBM = lens.apertureIndex (index)
%

%
%INPUT
%   obj: lens object of SCENE3D
%   index: specify the number of the surface which works as 'diaphrgam'
%   
%   
%
%OUTPUT
%
%   varargout: none for the moment
%
% MP Vistasoft 2014


%% CHECK FOR NUMERIC INPUT

if not(exist('index')) || not(isnumeric(index))
    error(' Specify a numeric value for INDEX')
end

%% GET OLD SURFACE
surf=obj.get('surfacearray',index);

%% SET the SUBTYPE as 'diaphragm'
surf.set('subtype','diaphragm')

%% Append to lens structure
obj.set('surface array index',surf,index)

