% Estimate the Phase Function of the pupil function for the specify coefficients and type of
% aberration

function [PhaseW]=psfGetPhaseFun(Coeff,type,ro,theta,varargin)

% INPUT
% Coeff: coeffs for the aberration
% type: type of aberration {'primary aberration','defocus'}
% ro: vector or matrix for radial coordinate
% theta: vector of matrix for angle coordinate
% varargin  {1}: specify if plot the result for DEBUG

% OUTPUT
% W: Phase function


%% Create a  MASK to limit the wavefront error within the pupil coordinate
Mask=circMask(ro);

%% Compute according to the type of coeffs
switch type
    case {'primary aberration';'pa';'seidel';'power series'}          
        fname=paGetFieldNames(Coeff);
        numF=length(fname); %number of coeffs
        
        %CHECK valued coeffs
        fLetter=(fname{1}(1));
        if length(fname{1})>3
            error ('Convert coeffs to include height field in the coeff value' )
        else
            switch fLetter
                case {'Z';'z'}
                    error (' Zernike coeffs NOT AVAILABLE yet!')
                case {'S';'s'}
                    error('Convert Seidel coeffs to Peak Coeffs')
            end
        end
        
        for fi=1:numF
            kC=str2num(fname{fi}(2)); %get radius order
            lC=str2num(fname{fi}(3)); %get cosine-theta order

            Wkl(:,fi)=getfield(Coeff,fname{fi});%weight for normalized coordinate            
            B(:,:,fi)=pa4thWaveAber(kC,lC,ro,theta); %bases
            % Number of wavelength
            nW=size(Wkl,1);
            for li=1:nW 
%                 WB(:,:,fi,li)=Wkl(li,fi).*B(:,:,fi);                 
                WB(:,:,fi,li)=Wkl(li,fi).*B(:,:,fi).*Mask;
            end
        end
        PhaseW=squeeze(sum(WB,3));

    case {'def';'defocus'}
        [fname]=paGetFieldNames(Coeff);
        numF=length(fname); %number of coeffs
        nW=size(getfield(Coeff,'wave'),1); %wavelength dependence 
%         nW=size(Coeff,1); %wavelength dependence 
        for ci=1:size(Coeff,2)
            value(:,ci)=getfield(Coeff,fname{ci});%weight for normalized coordinate
            B(:,:,ci)=pa4thWaveAber(ci*2,0,ro,theta); %bases
            for li=1:nW
              WB(:,:,ci,li)=value(li,ci).*B(:,:,ci).*Mask;
%                 WB(:,:,ci,li)=Coeff(li,ci).*B(:,:,ci).*Mask;
            end            
        end 
        PhaseW=squeeze(sum(WB,3));
    case {'defocusBraat';'defocus Braat';'def Braat'}
        [fname]=paGetFieldNames(Coeff);
        numF=length(fname); %number of coeffs
        nW=size(getfield(Coeff,'wave'),1); %wavelength dependence 
        for ci=1:size(Coeff,2)
            value(:,ci)=getfield(Coeff,fname{ci});%weight for normalized coordinate
            B(:,:,ci)=pa4thWaveAber((ci-1)*2,0,ro,theta); %bases
            for li=1:nW
%                 WB(:,:,ci,li)=Coeff(li,ci).*B(:,:,ci).*Mask;
                WB(:,:,ci,li)=value(li,ci).*B(:,:,ci).*Mask;
            end            
        end 
        PhaseW=squeeze(sum(WB,3));
        
    otherwise
        error(['Not valid ',type,' as type of aberration']) 
end

%% DEBUG
if nargin>4
    flagDebug=varargin{1};
else
    flagDebug=[];
end

%Plot

if not(isempty(flagDebug))&&( strcmp(flagDebug,'true')  || strcmp(flagDebug,'debug'))
    x_p=ro.*cos(theta); % x-axis normalized pupil coordinate
    y_p= ro.*sin(theta); % y-axis normalized pupil coordinate
    if nargin>5
        nW0=varargin{2}; %which wavelength
    else
        prompt=['Select a wavelength (specify the index): '];
        nW0=input(prompt);
    end
    surf(x_p,y_p,PhaseW(:,:, nW0))
%     contour(x_p,y_p,PhaseW(:,:, nW0))
%     contourf(x_p,y_p,PhaseW(:,:, nW0))
    xlabel('normalized x-axis '),ylabel('normalized y-axis'),zlabel('Wavefront')
    title(['Phase function for ', type, ' coefficients'])
    
end