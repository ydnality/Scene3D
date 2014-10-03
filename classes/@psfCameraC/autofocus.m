function varargout=autofocus(obj,wave0, waveUnit,varargin)
% Set the film of the camera in focus for the specified wavelength
%
%   psfCamera.autofocus(wave0, waveUnit)
%Examples
%   1:    psfCamera.autofocus(555, 'nm')
%   2:    psfCamera.autofocus(0.555, 'mm')
%   3:    psfCamera.autofocus(2, 'index') second value in the vector
%
% The camera has a point source, lens, and film.
%
% INPUT
% wave0
% waveUnit = 'nm' or 'mm'or 'index'
% varargin  {1}: n_ob
% varargin  {2}: n_im


% OUTPUT
% varargout{1}: new film position
%
%

% MP Vistasoft Team, Copyright 2014


%% GET wavelength vector
wave=obj.get('wave');

%% HAVE I t

if not(exist('wave0')) || not(exist('waveUnit'))
    error ('Specify the wavelength for the autofocus, example:  psfCamera.autofocus(555, "nm")')
end

switch waveUnit
    case {'nm'}
        wave0=wave0; %wave is in nm as well as wave0
    case {'um'}
        wave0=wave0*1e3; %wave is in nm, instead wave0 wan in um
    case {'mm'}
        wave0=wave0*1e6; %wave is in nm, instead wave0 wan in mm
    case {'m'}
        wave0=wave0*1e9; %wave is in nm, instead wave0 wan in m
    case {'index';'ind'}
        wave0=wave(wave0); % get the wave specified by the index
    otherwise
        error ('Specify a wavelength for the autofocus, example:  psfCamera.autofocus(555, "nm") ')
end

%index of the selected wavelength
ind0=find(wave==wave0);

%check if the selceted wavelenght exist and is just one 
if (isempty(ind0)) || (length(ind0)>1)
    error (' Not found a "unique" wavelength matching to the selected ones!!!!')
end


%% SET the film position in focus 
% Get input
% lens=obj.lens;
lens=obj.get('lens');
% film=obj.film;
film=obj.get('film');
% pSource=obj.pointSource;
pSource=obj.get('pointSource');

if nargin>3
   n_ob = varargin{1};    n_im = varargin{2};   
else
    n_ob = 1;    n_im = 1;                    
end

%Find Gaussian Image Point (wavelength dependence)
imagePoint = lens.findImagePoint(pSource,n_ob,n_im);

%Right distance for the selected wavelength
dist0=imagePoint(ind0,3); % z position


%Get previous film position
oldPos=obj.get('film');

% Set new distance
newPos=oldPos;
newPos.position(3)=dist0;

% film=pbrtFilmC('position', oldPos, 'size', filmSize, 'wave', wave, 'resolution', resolution);
%% SET THE NEW FILM
obj.set('film',newPos);

if nargout>0
    varargout{1}=dist0;
end
