% DRAW the source the specify point


function [out]=drawPoint(pointSource,wave0,wave,coord_type)

% INPUT
% Pupil: struct   .zpos  (optical axis position)
                      % .diam   (diameter)
% wave0: specify the wavelength to plot
% wave: set of all possible wavelength
% coord_type: specify which coordinate to plot {'x';'y'}
                      

%OUTPUT
%out: 0 or 1


%% CHECK if wavelength matches

if isempty(wave0) || isempty(wave)
    indW=1;
else
    wave0=550; %nm   select a wavelengt
    indW=find(wave==wave0);

    if isempty(indW)
        error(['Not valid matching between wavelength ',wave0,' among ',wave])
    end
end


%% Which coordinate

pZ=pointSource(indW,3); %Z coord
switch coord_type
    case {'x';'x-axis'}
        pH=pointSource(indW,1); 
    case {'y';'y-axis'}
        pH=pointSource(indW,2); 
    case {'eccentricity';'x-y';'r'}
        pH=sqrt(pointSource(indW,1).^2+pointSource(indW,2).^2);
        pAngle=atan(pointSource(indW,2)/pointSource(indW,1)) ; % useful for 3D plot
    otherwise
        error(['Not valid ',coord_type,' as axis'])
end


%% Parameters for the PLOT
colorLine='--g'; %black line
Lwidth=2; %line width
Msize=6; %mark size
Mcolor=''; %mark colour
%  Point source
hold on
stem(pZ,pH,colorLine,...
    'MarkerFaceColor','red','MarkerEdgeColor','blue','MarkerSize',Msize,'LineWidth',Lwidth)


%% SET OUTPUT
out=1;