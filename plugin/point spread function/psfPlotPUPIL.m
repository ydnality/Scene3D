% Graphic the  Wavefront Error (aka Pupil function) for the selected
% wavelength in 'real' coordinate

function [out]=psfPlotPUPIL(Pupil,xn,yn,ExP_Diam,limit,method,wave,wave0,varargin)

%INPUT
%
%xn:normalized  pupil coordinate
%yn: normalized  pupil coordinate
% ExP_diam: exit pupil diameter
%limit=[x_p limit, y_p limit]
%method: for plotting
%wave:
%wave0: selected wavelength

%OUTPUT

%% CHECK INPUT
inW0=find(wave==wave0);

if isempty(inW0)
    error('The select wavelength does NOT match with any sample')
else
    inW0=inW0(1);
end
%then 
x_p=xn.*ExP_Diam(inW0,:);y_p=yn.*ExP_Diam(inW0,:);

if isempty (limit)
    indx=1:length(x_p);
    indy=1:length(y_p);
elseif length(limit)==1
    indx=(abs(x_p)<=limit);
    indy=indx;    
else
    indx=(abs(x_p)<=limit(1));
    indy=(abs(y_p)<=limit(2));
end

%% CREATE MESH FOR PLOT
[X,Y]=meshgrid(x_p(indx),y_p(indy));
waveFront=angle(Pupil(indx,indy,inW0));

%% PLOT
%TITLE
FontSize1=12;

if nargin<8
    figure
end

switch method
    case {'surf'},surf(X,Y,waveFront)
    case {'contour'},contour (X,Y,waveFront)
    case {'contour3'},contour3 (X,Y,waveFront)
    case {'contourf'},contourf (X,Y,waveFront)
       
end
% set(gca,'LineWidth',LineWidth,'LineStyle','--')
% colormap hot
colorbar
title(['Exit Pupil Wavefront (Phase) Error at ',num2str(wave0*1e6),'nm'])
xlabel(['x [mm]']),ylabel(['y [mm]'])
zlabel('# of wavelength /rad')

out=1;


