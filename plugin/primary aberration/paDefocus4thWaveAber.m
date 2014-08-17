%Compute the Defocus Coeffs for the normalized coordinate
% based on different level of approxiamtion of Nijboer-Zernike Theory
% ADVICE: select 'highNA' as method to compute NA

function [W]=paDefocus4thWaveAber(zI,zG,NA,ro,varargin)

%INPUT
%zI: distance of image plane
%l: distance of Gaussian image point
%NA: numerical aperture
%ro:  normalized radius coordinate
%vararging {1}: method of approximation


%OUTPUT
%W: Term of Defocused 4th order Wavefront aberration

%NOTE:
% NA, zI and zG has to be scalar not allow to be wavelength dependent

%% See APPENDIX A: 
%%Braat, Joseph, Peter Dirksen, and Augustus JEM Janssen. "Assessment of an extended Nijboer–Zernike approach for the computation of optical point-spread functions." JOSA A 19.5 (2002): 858-870.
% Wf=dZ(1-sqrt(1-sin(p)^2); p=ro*NA


if nargin>4
    approx_method=varargin{1};
else
    approx_method='smallNA'; %small Numerical Aperture for default
end

%% DEFOCUS ALONG Z axis
dZ=zG-zI; %defocus Z
%%
switch approx_method
    
    case {'smallNA';'small NA';'small'}
        %% FOR SMALL NUMERICAL APERTURE        
        W=0.5*dZ.*NA.^2*ro.^2;% output
    case {'highNA';'high NA';'high'}
        %% FOR HIGH NUMERICAL APERTURE
        % 1-sqrt(1-sin(p).^2) approx as quadratic poly b0+b1*ro^2
        %for small NA b0=0; b1=
        sinP= NA.*ro;
        val=1-sqrt(1-sinP.^2);
        %Fitting with second order poly
        p=polyfit(ro.^2,val,2);
        %Evaluate such poly
        f=polyval(p,ro.^2);
        W=dZ.*f;%output
        
    case {'debug'}
        %% FOR SMALL NUMERICAL APERTURE        
        Wsmall=0.5*dZ.*NA.^2*ro.^2;% output
         %% FOR HIGH NUMERICAL APERTURE
        % 1-sqrt(1-sin(p).^2) approx as quadratic poly b0+b1*ro^2
        %for small NA b0=0; b1=
        sinP= NA.*ro;
        val=1-sqrt(1-sinP.^2);
        %Fitting with second order poly
        p=polyfit(ro.^2,val,2);
        %Evaluate such poly
        f=polyval(p,ro.^2);
        Whigh=dZ.*f;%output
        %Difference
        Err=Whigh-Wsmall;
        Err_rel=Err./Whigh*100; %percentual error
        %elimate unchanged value due to polyfit
        ind=find(Err_rel==100);
        Err_rel(ind)=0;
       %Plot 
       figure
       subplot(2,2,1)
       surf(Wsmall)
       title('Wavefront Aberration: Small Numerical Aperture')
       subplot(2,2,2)
       surf(Whigh)
       title('Wavefront Aberration: High Numerical Aperture')
       subplot(2,2,[3 4])
       surf(Err_rel)
       title('Wavefront Aberr. Relative Error between HighNA and SmallNA')
       zlabel('% error ')
       %OUTPUT
       W=Whigh;
        
    otherwise
        error (['Not accepted ',approx_method,' as method to approx defocus!'])
end
