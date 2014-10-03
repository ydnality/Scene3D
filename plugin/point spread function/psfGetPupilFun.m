% Build Pupil function according to apodization, phase function and
% wavelegth


function [PupilFun]=psfGetPupilFun(ApodW,W_PA,W_Def,wave,varargin)


%INPUT
% ApodW: apodization function
% W_PA: Phase function for primary aberration
% W_Def: Phase function for defocus
% wave: vector wavelength
% varargin  {1}: flag for debug  {2} debug parameter (wavelength to plot)

% OUTPUT
%PupilFun: wavelength-dependent pupil function



%% CHECK if all input are wavelength dependent

nW=size(wave,1); %wavelength number

if not(size(ApodW,3)==nW)
    ApodW=repmat(ApodW,[1,1,nW]);
end
if not(size(W_PA,3)==nW)
    W_PA=repmat(W_PA,[1,1,nW]);
end
if not(size(W_Def,3)==nW)
    W_Def=repmat(W_Def,[1,1,nW]);
end



%% Build the pupil function for each wavelenght

for li=1:nW
    K=-i*2*pi/wave(li,1); % pahse converison factor  2*pi/lambda
    PhaseFun=K.*(W_PA(:,:,li)+W_Def(:,:,li)); %Phase function
    PupilFun(:,:,li)=ApodW(:,:,li).*exp(PhaseFun); %Pupil function
end

%% DEBUG

if nargin>4
    switch varargin{1}
        case {'debug'}
            if nargin>5
                wave0=varargin{2};
                indW0=find(wave==wave0);
            else
               prompt=['Select a wavelength (specify the index), among ',num2str(wave'*1e6),' [nm]: '];
                indW0=input(prompt);
                wave0=wave(indW0);
            end
            if nargin>6
                typeData=varargin{3};
            else
                typeData='abs';
            end
            switch typeData
                case {'modulus';'abs'}
                    surf(abs(PupilFun(:,:,indW0)))
                    title(['Pupil Function Modulus at ',num2str(wave0*1e6),' [nm]']) 
                case {'angle';'phase'}
                    surf(angle(PupilFun(:,:,indW0)))
                    title(['Pupil Function Phase at ',num2str(wave0*1e6),' [nm]'])
            end
                       
        otherwise
            warning('Not valid input for DEBUG')
    end
end
