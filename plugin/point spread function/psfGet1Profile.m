%Get a 1D profile of a 2D PSF according to the selected axis


function [PSF1D,varargout]=psfGet1Profile(PSF2D,xD,yD,wave,wave0,method,varargin)

%INPUT
%PSF2D: [MxL]matrix for a 2D point spread function
%xD: [NxM] 
%yD: [NxL] 
%wave:[Nx1]
%wave0:
%method: select method to obtain profile


%OUTPUT
% PSF1F: 1D profile
%varargin: coordinate for evaluation

inW0=find(wave==wave0);

if isempty(inW0)
    error('The select wavelength does NOT match with any sample')
else
    inW0=inW0(1);
end


switch method
    
    case {'peak-x'}
        %find peak
        peakV=max(max(PSF2D(:,:,inW0)));
        [rP,cP]=find(PSF2D(:,:,inW0)==peakV);
        if numel(rP)>1
            %In case of multiple peak select the median position
            in0=find(rP==ceil(median(rP)));
            rP=rP(in0(1));
            cP=cP(in0(1));
        end
        % check if available sampling limit
        if nargin>6
            x1d=varargin{1};
        else
            x1d=xD(inW0,:);
        end
        PSF1D=interp1(xD(inW0,:),PSF2D(rP,:,inW0),x1d);
        if nargout>1
            varargout{1}=x1d;
        end
        
    case {'peak-y'}
        %find peak
        peakV=max(max(PSF2D(:,:,inW0)));
        [rP,cP]=find(PSF2D(:,:,inW0)==peakV);
        % check if available sampling limit
        if nargin>6
            y1d=varargin{1};
        else
            y1d=yD(inW0,:);
        end
        PSF1D=interp1(yD(inW0,:),PSF2D(:,cP,inW0),y1d);
        if nargout>1
            varargout{1}=y1d;
        end
end

