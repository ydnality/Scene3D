% Graphic the PSF2D for the given wavelength


function [out]=psfPLOT(PSF,x_im,y_im,limit,method,wave,wave0,varargin)

%INPUT
%
%
%limit=[x_im limit, y_im limit]
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
x_im=x_im(inW0,:);y_im=y_im(inW0,:);

if isempty (limit)
    indx=1:length(x_im);
    indy=1:length(y_im);
elseif length(limit)==1
    indx=(abs(x_im)<=limit);
    indy=indx;    
else
    indx=(abs(x_im)<=limit(1));
    indy=(abs(y_im)<=limit(2));
end

%% CREATE MESH FOR PLOT
[X,Y]=meshgrid(x_im(indx),y_im(indy));
PSF=PSF(indy,indx,inW0);

%% PLOT
%TITLE
FontSize1=12;

if nargin<8
    figure
end

switch method
    case {'surf'},surf(X,Y,PSF)
    case {'contour'},contour (X,Y,PSF)
    case {'contour3'},contour3 (X,Y,PSF)
    case {'contourf'},contourf (X,Y,PSF)
       
end
% set(gca,'LineWidth',LineWidth,'LineStyle','--')
% colormap hot
colorbar
title(['PSF at ',num2str(wave0*1e6),'nm'])
xlabel(['x [mm]']),ylabel(['y [mm]'])
zlabel('Normalized intensity')

out=1;


