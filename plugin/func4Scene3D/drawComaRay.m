% DRAW the coma rays from the point source to the rims of the specify pupil


function [out]=drawComaRay(pointSource,Pupil,wave0,wave,coord_type)

% INPUT
% Pupil: struct   .zpos  (optical axis position)
                      % .diam   (diameter)
% wave0: specify the wavelength to plot
% wave: set of all possible wavelength
% coord_type: specify which coordinate to plot {'x';'y'}
                      

%OUTPUT
%out: 0 or 1




%% CHECK if wavelength matches

wave0=550; %nm   select a wavelengt
indW=find(wave==wave0);

if isempty(indW)
    error(['Not valid matching between wavelength ',wave0,' among ',wave])
end

%% Which coordinate

if  size(pointSource,1)>1
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
else
    pZ=pointSource(3); %Z coord
    switch coord_type
        case {'x';'x-axis'}
            pH=pointSource(1); 
        case {'y';'y-axis'}
            pH=pointSource(2); 
        case {'eccentricity';'x-y';'r'}
            pH=sqrt(pointSource(1).^2+pointSource(2).^2);
            pAngle=atan(pointSource(2)/pointSource(1)) ; % useful for 3D plot
        otherwise
            error(['Not valid ',coord_type,' as axis'])
    end
end



%% PUPIL
P_zpos=Pupil.zpos(indW); %Z coord
P_upH=Pupil.diam(indW)/2; % Upper rim of the pupil
P_loH=-Pupil.diam(indW)/2; % lower rim of the pupil


%% Parameters for the PLOT
colorLine='--g'; %black line
Lwidth=2; %line width

% Upper Coma Ray
hold on
plot([pZ P_zpos],[pH P_upH],colorLine,'LineWidth',Lwidth)
% LowerComa Ray
plot([pZ P_zpos],[pH P_loH],colorLine,'LineWidth',Lwidth)

%% SET OUTPUT
out=1;