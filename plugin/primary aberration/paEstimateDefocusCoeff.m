

function [Coeff,varargout]=paEstimateDefocusCoeff(zI,zG,NA,n_im,approx_method,varargin)
%Compute the Defocus Coeffs for the normalized coordinate
% based on different level of approxiamtion of Nijboer-Zernike Theory
% ADVICE: select 'highNA' as method to compute NA
%
% function [Coeff]=paEstimateDefocusCoeff(zI,zG,NA,n_im,approx_method,varargin)
%
%INPUT
%zI: position of image plane
%zG: position of Gaussian point
%l: distance of Gaussian image point
%NA: numerical aperture
%n_im: refractive index in the image space
%vararging {1}: method of approximation

%OUTPUT
%W: Term of Defocused 4th order Wavefront aberration
% varargout {1}: defocus distance

%NOTE:
% NA, zI and zG has to be scalar not allow to be wavelength dependent
%
% MP Vistasoft 2014

%% See APPENDIX A: 
%%Braat, Joseph, Peter Dirksen, and Augustus JEM Janssen. "Assessment of an extended Nijboer–Zernike approach for the computation of optical point-spread functions." JOSA A 19.5 (2002): 858-870.
% Wf=dZ(1-sqrt(1-sin(p)^2); p=ro*NA
% 
%and 'The Extended Nijboer-Zernike Driffraction Theory and its
%Applications' Sven van Haver Thesis


if nargin>5
    ro=varargin{1};
else
    % default parameters
    nSample=2^6;
    ntime=1;
    
    range_pupil=2*ntime ; % al range od pupil function in normalized coordinate
    d_pupil=range_pupil/nSample; %sampli
    xn=[-nSample/2:nSample/2-1]*d_pupil;
    yn=[-nSample/2:nSample/2-1]*d_pupil;
    [Xn,Yn]=meshgrid(xn,yn);
    ro=circMask(sqrt(Xn.*Xn+Yn.*Yn)).*sqrt(Xn.*Xn+Yn.*Yn); 
%     ro=(sqrt(Xn.*Xn+Yn.*Yn));  
end

%% DEFOCUS ALONG Z axis
dZ=zG-zI; %defocus Z
%%
switch approx_method
    
    case {'smallNA';'small NA';'small';'1term'}
        %% FOR SMALL NUMERICAL APERTURE  
        Coeff1=0.5.*(NA./n_im).^2.*dZ; %Unique coeff
%         W=0.5*dZ.*NA.^2*ro.^2;% output
        Coeff=Coeff1;
    case {'2term';'2terms'}
        s0=NA./n_im; % See definition vanHaver's thesis- pag 22
        C1=(s0.^2)/2; % coeff for  ro^2
        C2=(s0.^4)/8; % coeff for ro^4
        Coeff=dZ*[C1 C2]; %coeffs with defocus
    case {'3term'; '3terms'}
        s0=NA./n_im; % See definition vanHaver's thesis- pag 22
        coeff2t=paEstimateDefocusCoeff(zI,zG,NA,n_im,'2term',ro); %coeffs up to second terms
        C3=(s0.^6)/16; % coeff for  ro^6
        Coeff=[coeff2t C3.*dZ]; %coeffs with defocus
    case {'highNA';'high NA';'high'}
        %Method proposed by Braat, Joseph, Peter Dirksen, and Augustus JEM Janssen. "Assessment of an extended Nijboer–Zernike approach for the computation of optical point-spread functions." JOSA A 19.5 (2002): 858-870.
        % Wf=dZ(1-sqrt(1-sin(p)^2); p=ro*NA
        %% FOR HIGH NUMERICAL APERTURE
        % 1-sqrt(1-sin(p).^2) approx as quadratic poly b0+b1*ro^2
        %for small NA b0=0; b1=
        sinP= NA./n_im.*ro;
        val=1-sqrt(1-sinP.^2);
        %Fitting with second order poly
        [p,S]=polyfit(ro.^2,val,1);
        %Evaluate such poly
%         f=polyval(p,ro.^2);
%         W=dZ.*f;%output
        Coeff=dZ.*p;
    case {'best'}
        %Expected value
        sinP= NA./n_im.*ro;
        val=1-sqrt(1-sinP.^2);
        % Methods
        list={'1term';'2term';'3term'};
        for ti=1:length(list)
            [C{ti}]=paEstimateDefocusCoeff(zI,zG,NA,n_im,list{ti},ro);
            valC{ti}=zeros(size(ro,1),size(ro,2));
            for ci=1:length(C{ti})
                exp=str2num(list{ti}(1))*2;
                valC{ti}=valC{ti}+C{ti}(ci).*ro.^exp;
            end
            % Residual
            res(ti)=sum(sum((val-valC{ti}).^2));
        end
        if NA<0.6 %threshold given by vanHaver's thesis- pag 22
            Coeff=C{1}; % only first term
        elseif (res(1)==0)
            Coeff=C{1}; % first term is enought
        else
            % Use more terms
            res0=res(1);
            perRes=(res(2:end)-res0)./res0; % relative residual 
        end
    case {'debug'}
        %% FOR SMALL NUMERICAL APERTURE 
        Coeff1=0.5.*NA./n_im.^2; %Unique coeff
        Wsmall=Coeff1.*dZ.*ro.^2;% output
         %% FOR HIGH NUMERICAL APERTURE
        % 1-sqrt(1-sin(p).^2) approx as quadratic poly b0+b1*ro^2
        %for small NA b0=0; b1=
        sinP= NA./n_im.*ro;
        val=1-sqrt(1-sinP.^2);
        %Fitting with second order poly
        [p,S]=polyfit(ro.^2,val,1);
        %Evaluate such poly
        Coeff2=p.*dZ;
        f=polyval(Coeff2,ro.^2);
        Whigh=f;%output
        %Difference
        Err=Whigh-Wsmall;
        % Err relative
        [indR,indC]=find(Wsmall~=0);
        Err_rel=zeros(size(Wsmall,1),size(Wsmall,2));
        T1=Err(indR,indC);
        T2=Wsmall(indR,indC);
        Err_rel(indR,indC)=T1./T2*100;
%         Err_rel=Err./Whigh; %percentual error
        %elimate unchanged value due to polyfit
        indINF=isinf(Err_rel);
        Err_rel(indINF)=0;
       %Plot 
       figure
       subplot(2,2,1)
       imagesc(Wsmall)
       title(['Wavefront Aberration: Small NA, ',num2str(NA)])
       colorbar
       subplot(2,2,2)
       imagesc(Whigh)
       title(['Wavefront Aberration: High NA, ',num2str(NA)])
       colorbar
       subplot(2,2,3)
       imagesc(Err)
       title('Wavefront Aberr.  Error between HighNA and SmallNA')
       colorbar
       zlabel('% error ')
       subplot(2,2,4)
       imagesc(Err_rel)
       title(' Relative Error (%) between HighNA and SmallNA')
       colorbar
       zlabel('% error ')
       %OUTPUT
       Coeff=p.*dZ;
        
    otherwise
        error (['Not accepted ',approx_method,' as method to approx defocus!'])
end


%% SET OUTPUT
if nargout>1
    varargout{1}=dZ;
end