% Estimate the PSF for an Imaging System according to its object and the
% film

function [PSF,x_im,y_im,varargout]=psfEstimate(ImagSyst,nSample,ntime)


%INPUT: to be completed
%ImagSyst:
%nSample: number of sample in each direction [power of 2]
%ntime: semi-range of definition for normalized pupil function (radius=1)


%OUTPUT: to be completed
%PSF
%x_im: vector of coordinate in image plane [wavelength dependent]
%y_im: vector of coordinate in image plane [wavelength dependent]


%% Get relevant parameters:
wave=ImagSyst.wave;
nW=size(wave,1);
n_im=ImagSyst.n_im;

% Defocus given by the difference between film position and ideal focus
if size(ImagSyst.film{end}.z_pos,1)==nW
    defocusZ(:,1)=ImagSyst.film{end}.z_pos; %last film
else
    defocusZ(:,1)=repmat(ImagSyst.film{end}.z_pos,nW,1);
end
if isinf(ImagSyst.object{end}.z_pos)
    defocusZ(:,2)=ImagSyst.cardPoints.dFi+ImagSyst.cardPoints.lastVertex; 
else
    defocusZ(:,2)=ImagSyst.object{end}.ConjGauss.z_im;
end

%Peak Value ABERRATION COEFFs
PeakCoeff=ImagSyst.object{end}.Wavefront.PeakCoeff;


%% GENERATE PUPIL FUNCTION
%Normalized coordinate of the pupil function (aperture goes form -1 to 1)
range_pupil=2*ntime ; % al range od pupil function in normalized coordinate
d_pupil=range_pupil/nSample; %sampling
%
xn=[-nSample/2:nSample/2-1]*d_pupil;
yn=[-nSample/2:nSample/2-1]*d_pupil;
% x=[-Rmax:dx:Rmax];y=[-Rmax:dy:Rmax];
% [Xn,Yn]=meshgrid(fx,fy);
[Xn,Yn]=meshgrid(xn,yn);

%Image plane  normalized coordinate 
d_image=1/range_pupil;
image_range=1/d_pupil;

xn_im=[-nSample/2:nSample/2-1]*d_image;
yn_im=[-nSample/2:nSample/2-1]*d_image;


%% COMPUTE PUPIL FUNCTION AND POINT SPREAD FUNCTION

% CREATE MASK
%Assumed uniform
Mask=circMask(sqrt(Xn.*Xn+Yn.*Yn));

% CREATE APODIZATION FUNCTION
Apod=circMask(sqrt(Xn.*Xn+Yn.*Yn));
    

%% CREATE PHASE WAVEFRONT ABERRATION
fname=paGetFieldNames(PeakCoeff);
numF=length(fname); %number of coeffs
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

switch wvf_type
    case {'zernike'}
        
    case {'peak'}
        ro=Mask.*sqrt(Xn.*Xn+Yn.*Yn); %normalized radius
        theta=Mask.*atan2(Yn,Xn); % theta  
%         theta=atan(Yn./Xn); % theta  
    
        for fi=1:numF
            kC=str2num(fname{fi}(2)); %get radius order
            lC=str2num(fname{fi}(3)); %get cosine-theta order
            
            Wkl(:,fi)=getfield(PeakCoeff,fname{fi});%weight for normalized coordinate            
            B(:,:,fi)=pa4thWaveAber(kC,lC,ro,theta); %bases            

            for li=1:nW 
                WB(:,:,fi,li)=Wkl(li,fi).*B(:,:,fi);  
            end
        end
        PhaseW=squeeze(sum(WB,3));
        %
end
        
        

%% DEFOCUS TERM
% used defocused pupil function described by Mahajan- PartII Wave
% Diffraction Optics, equation 1-47



for li=1:nW
    ExPDiam(li,1)=ImagSyst.object{end}.Radiance.ExP.diam(li,1)-ImagSyst.object{end}.Radiance.ExP.diam(li,2);
%     efl(li,1)=ImagSyst.cardPoints.fi(li,1);
    efl(li,1)=ImagSyst.object{end}.ConjGauss.z_im(li,1)-mean(ImagSyst.object{end}.Radiance.ExP.z_pos(li,:),2);
    [NA(li,:)]=paraxNumAperture(ExPDiam(li,1),efl(li,1),n_im(li,1));
    PhaseDefocus(:,:,li)=paDefocus4thWaveAber(defocusZ(li,1),defocusZ(li,2),NA (li,:),ro,'small');
%     PhaseDefocus(:,:,li)=paDefocus4thWaveAber(defocusZ(li,1),defocusZ(li,2),NA (li,:),ro,'debug');
end


%% COMPUTE PSF FOR ALL TERM

%Phase term
Kwave=-1i*2*pi./wave;

for li=1:nW
    Pupil(:,:,li)=Apod.*exp(Kwave(li,1).*(PhaseW(:,:,li)+PhaseDefocus(:,:,li)));
    PSF(:,:,li)=psfPupil2PSF(Pupil(:,:,li),'incoherent');
    
    %and its coordinate
    
    [x_im(li,:),y_im(li,:)]=psfNormalized2RealCoordinate(xn_im,yn_im,wave(li,1),NA(li));
end



%% SET OUTPUT

PSF=PSF;
x_im=x_im;
y_im=y_im;

if nargout>3
    varargout{1}=Pupil; %Pupil complex function
    if nargout>4
        varargout{2}=ExPDiam; %ExP diamter
        if nargout>5
            varargout{3}=NA; %numerical aperture
        end
    end
end
