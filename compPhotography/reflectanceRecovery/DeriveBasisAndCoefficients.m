% DeriveBasisAndCoefficients.m
%
% Derive the coefficients from a multispectral capture.  The original data
% are a matrix with, say, 6 color channels based on the D100 with a clear
% image and with a red filter.  We compute color signal basis functions (up to 6).
% We then use the camera data to estimate the values of the basis
% coefficients.  The coefficients and basis functions are stored in a file
% that can be called by vCamera and other programs so that we have
% estimated color signals.
%
% Feng Xiao, Brian Wandell and Joyce Farrell
%
% mc stands for multiple capture
% ms stands for multi-spectral
% lm stands for linear model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Step 1: Select basis functions
wavelength = 400:10:700;
% lightList = {'A','B','C','D50','D55','D65','D75','FL11','FL2','FL7',...
%         'SimonFraserIlluminants','OfficeFL','Vivitar'};
lightList = {'A','B','C','D50','D55','D65','D75'};

lights = lmLookupSignals(lightList, 'illuminants', wavelength,1);
% plot(lights)

surfaceList = {'Clothes','Food','Hair','Objects','Nature','Paint','macbethChart'};
% This one has a NaN in it.  surfaceList = {'SkinReflectance'}
surfaces = lmLookupSignals(surfaceList, 'surfaces', wavelength,0);

nBases = 4;
[basis,cs,sValues] = lmColorSignalBasis(lights,surfaces,nBases);
csBasis.basis = basis;
csBasis.wave = wavelength;

plot(csBasis.wave,csBasis.basis,'-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Step 2: Compute the basis coefficients using the camera data
%  The data files read in here are created using the script
%  mcCombineExposureColor
%
% Read in a saved data file from that script.

partialName = fullfile('April22','DSC_0002-DSC_0005.mat')
load(partialName)

% We need to convert this to vcReadSpectra() format.
D100 = specQuery('sensors','D100',wavelength)/10000;
filters = specQuery('62mmFilters','Tiffen Red 29',wavelength);
sensor = [D100 D100.*repmat(filters,[1 3])];

msCoef = mcCamera2CSBasis(sensor, csBasis.basis, camRGB);  % this line computes the coefficients

partialName = fullfile('April22','macbeth1.mat')
mcSaveCoefAndBasis(partialName,msCoef,csBasis,vcInfo)

% To check the spd, try this:
% spd = rgbLinearTransform(msCoef,csBasis');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Step 5: display the constructed multispectral image on a sRGB monitor 
%  phosphor = getSensorSpectral('SRGB',wavelength);
%  I think we should probably be using the xyz2srgb calls in vCamera or IE.
%   For now, though, we are doing it Feng's way.
%
% XYZspectral = specQuery('sensors','CIEXYZ',wavelength);
% phosphors = XYZspectral * [ 3.2406   -0.9689    0.0557;  -1.5372  1.8758   -0.2040;   -0.4986    0.0415    1.0570];
% plot(phosphors)

gam = 2; RGB = mcDisplay(msCoef,csBasis,[],gam);
imagescRGB(RGB);

% Work stopped here ....

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Step 6: find the error between the estimation and the true spectral data
%  of macbeth color patches

% position of macbeth color patches
[x,y] = meshgrid(round(linspace(648,1076,6)),round(linspace(493,734,4))); 
x=x'; x=x(:); y=y'; y=y(:);

for i=1:length(x)
    patch = coef(:,y(i)-3:y(i)+3,x(i)-3:x(i)+3); 
    ef(:,i) = squeeze(mean(mean(patch,3),2));
end

est = colorSignalBasis*ef;

% load measured spectral data
load('macbethFL_TG.mat');
macbethTG
figure; plot(wavelength,est);
i=1;
