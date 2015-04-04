% Computer 1D PSF  along the directional-axis of a 2D Light Field

function [PSF2D]=LFComputePSF2D(PSF1D,method,varargin)


%INPUT
%PSF1D: psf along 1 axis (contains .Y,.normI)
%method: specify the method of computation for the 2D PSF
%          -> 'meridional' varargin{1} 
error('Deprecated');

%OUTPUT
%PSF2D: is the point spread funtion along the axis of the given LF object 
%       

switch method
    case {'meridional';'standard'}
        %Create OUTPUT sampling space
        v_Y=PSF1D.Y;
        [X,Y]=meshgrid(v_Y,v_Y); %create mesh
        R=(X(:).^2+Y(:).^2).^.5;
        for li=1:size(PSF1D.normI,1)
            %Append relevant data
            v_normI=PSF1D.normI(li,:);        
            PSF2d=X;
            PSF2d(:)=interp1(v_Y,v_normI,R); 
            % OUTPUT
            MnormI=reshape(PSF2d,size(Y,1),size(Y,2));
            %% CORRENCT NaN solution
            [indNaN]=find(isnan(MnormI));
            MnormI(indNaN)=0;
            OUTnormI(:,:,li)=MnormI;
        end
        %Set output
        PSF2D.X=X;PSF2D.Y=Y;
        PSF2D.normI=OUTnormI;
        
    case {'sagittal'}
        warning (['Method ',method,' has to be completed!!'])
        PSF2D=[];
    otherwise
        warning ('Not valide type of computation for PSF2D')
        PSF2D=[];
end





