%% COMPUTE THE PHASE WAVEFRONT FUNCTION for the given coeffs



function [PhaseW,ro,varagout]=paComputePhaseWvf(PeakCoeff,xn,yn,Mask)

%INPUT



%OUTPUT


%% INITIAL PARAMETER


[Xn,Yn]=meshgrid(xn,yn);

%STEP 1: check type of available coeffs
fname=paGetFieldNames(PeakCoeff);
numF=length(fname); %number of coeffs
nW=size(getfield(PeakCoeff,fname{1}),1);
%check is available Zernike or Peak Coeffs
switch fname{1}(1)
    case {'C','c'}
        wvf_type='zernike';
    case {'W','w','A','a'}
        if length(fname{1})==3
            wvf_type='peak';
        else
            error([fname{1},' is not a valid type of aberration. Probabily the peak value has to be computed '])
        end
    otherwise
        error([fname{1},' is not a valid type of aberration '])
end


%GENERATA WAVEFRONT
switch wvf_type
    case {'zernike'}
        
    case {'peak'}
        ro=Mask.*sqrt(Xn.*Xn+Yn.*Yn); %normalized radius
        theta=Mask.*atan2(Yn,Xn); % theta  
    
        for fi=1:numF
            kC=str2num(fname{fi}(2)); %get radius order
            lC=str2num(fname{fi}(3)); %get cosine-theta order
            
            Wkl(:,fi)=getfield(PeakCoeff,fname{fi});%weight for normalized coordinate           
            B(:,:,fi)=pa4thWaveAber(kC,lC,ro,theta); %bases            

            for li=1:nW
                WB(:,:,fi,li)=Wkl(li,fi).*B(:,:,fi);  
%                 WB(:,:,fi,li)=Wkl_real(li,fi).*B(:,:,fi);  
            end
%             %DEBUG
%             inW=3; %3th wavelength sample
%             wave0=wave(3);% select wavelength 55nm
        end
end
        PhaseW=squeeze(sum(WB,3));
%% OUTPUT
PhaseW=PhaseW;
ro=ro;
%OTHER OUTPUT
if nargout>2
    varagout{1}=WB;
    if nargout>3
        varargout{1}=B;
    end
end