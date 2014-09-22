% Graphic the PSF2D for the given wavelength


function [out]=psfPLOT(PSF,x_im,y_im,limit,method,wave,wave0,varargin)

%INPUT
%PSF : wavelength dependent PSF [Mx x My x N]
%x_im: wavelength dependent x coordinate vector [N x Mx] or matrix (meshed)[My x Mx x N]
%y_im: wavelength dependent y coordinate vector [N x My] or matrix (meshed)[My x Mx x N]
%limit=[x_im limit, y_im limit]
%method: specify the plot graph type
%wave: number of wavelength [N]
%wave0: selected wavelength

%OUTPUT
%

%% CHECK INPUT
inW0=find(wave==wave0);

if isempty(inW0)
    error('The select wavelength does NOT match with any sample')
else
    inW0=inW0(1);
end
%then CHECK coordinate dimensions


if not(ndims(x_im)==ndims(y_im))
    error (' INPUT dimensions for x- and y- coordinate DO NOT MATCH !!')
elseif ndims(x_im)==2 %y_im dimension is the same of x_im
    % input coordinate dimension  x_im=[N,Mx]; y_im=[N,My]
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
    
elseif ndims(x_im)==3 %y_im dimension is the same of x_im
    % input coordinate dimension  x_im=[My,Mx,N]; y_im=[My,Mx,N]
    x_im=x_im(:,:,inW0);y_im=y_im(:,:,inW0);
    %Check homogenity in input coordinate dimension
    if not(ndims(x_im)==ndims(y_im))
        error(['input coord dimension DO NOT MATCH, x_im: ',num2str(ndims(x_im)),' y_im: ',num2str(ndims(y_im))])
    end
    
   
    if isempty (limit)
        indx=1:size(x_im,2);
        indy=1:size(y_im,1);
    elseif length(limit)==1
        [indy0,indx]=find(abs(x_im)<=limit);
        [indy,indx0]=find(abs(x_im)<=limit);   
    else
        [indy0,indx]=find(abs(x_im)<=limit(1));
        [indy,indx0]=find(abs(x_im)<=limit(2)); 
    end
    indy=squeeze(indy); indx=squeeze(indx);
    %% CREATE MESH FOR PLOT
    X=x_im(indy,indx);Y=y_im(indy,indx);    
    
else
    error ('the dimension of the input coordinate exceed the limit (2 or 3) !')
end


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


