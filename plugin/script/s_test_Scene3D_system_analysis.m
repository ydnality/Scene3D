%s_test_Scene3D_system_analysis.m

%Test the conversion of a system from SCENE 3D to format of PSF3D


%% Initialize Iset
s_initISET


%% STEP 1: CREATE OPTICAL SYSTEM
%% Declare ray-trace type

% rtType = 'realistic';  %ideal/realistic
% debugLines = 50;
% % (lowerLeftx,LowerLefty,uppRightx,upperRighty)
% % Percentage of rectangular lens, 0 is middle
% subSection = [];   % Whole thing
% subSection = [-.25 -.25 .25 .25];   % Not working
%subSection = [0 0 .25 .25];   

    %% 1.1 Declare point sources
% declare point sources in world space.  The camera is usually at [0 0 0],
% and pointing towards -z.  We are using a right-handed coordinate system.
pointSource = [10 20 -1000];

%% Declare camera properties

% Build a sensor (film) object
% Position, size,  wave, waveConversion, resolution
% film = pbrtFilmC([0 0 51.2821	],[.2/sqrt(2) .2/sqrt(2)], 400:10:700, [(400:10:700)' (1:31)'], [50 50 31]);

% Declare film
filmPosition = [0 0 51.2821	];  % Good for 2Elens
% filmPosition = [0 0 37.4];        % Good for dgauss.50mm.  True focal about 37.3mm
% filmPosition = [0 0 107]; 

filmDiag = 3;  % Millimeters
filmSize = [filmDiag/sqrt(2) filmDiag/sqrt(2)];

wave = 400:50:600;
resolution =  [300 300 length(wave)];
film = pbrtFilmC('position', filmPosition, 'size', filmSize, 'wave', wave, 'resolution', resolution);

% Declare Lens
diffractionEnabled = false;   %diffraction causes imaginary directions!! TODO:  INVESTIGATE THIS!!
apertureDiameterMM = 2.2727;  %f/22 
% apertureDiameterMM = 3.1250;  %f/16
% apertureDiameterMM = 4.5455;  %f/11
fLength = 50;
apertureSamples = [51 51];
name = 'idealLensTest';
type = 'idealLens';
jitterFlag = true;

thinlens = lensC('name', name, 'type', type, 'focalLength', fLength, 'diffractionEnabled', diffractionEnabled, 'wave', wave, 'aperturesample', apertureSamples);

% Make a function that goes gets a lens from a file, say
% lens = lensC('filename',fname);
%  That should go and see if there is a file called fname.
%
lensFile = fullfile(s3dRootPath, 'data', 'lens', '2ElLens.mat');
import = load(lensFile,'lens');
thickLens = import.lens;
thickLens.apertureMiddleD = 4;

% lensFile = fullfile(s3dRootPath, 'data', 'lens', 'dgauss.50mm.mat');
% import = load(lensFile,'lens');
multiLens = import.lens;

% lens = multiLens;
lens = multiLens;
lens.set('wave', wave);
% Matrix of n for each surface element.
% Apertures are 0
n = lens.get('nArray');

%% CREATE IMAGING SYSTEM
ppsfCamera = ppsfCameraC('lens', lens, 'film', film, 'pointSource', pointSource);

%% CALL Function to conversion

[result1]=paraxAnalyzeScene3DSystem('lens',lens);

[result2]=paraxAnalyzeScene3DSystem('psfCamera',ppsfCamera);

[result3]=paraxAnalyzeScene3DSystem('specify',lens,film,pointSource);

%% FIND the image point

Hobj=result2.cardinalPoint.ObjectSpace.principalPoint; %principal point in the object space
Him=result2.cardinalPoint.ImageSpace.principalPoint; %principal point in the image space
n_ob=1; %refractive index in object space
n_im=1; %refractive index in image space
focalLength=result2.focallength; %effective focal length
magn=result2.magnification.lateral; %lateral magnification
abcdMatrix=result2.abcdMatrix; %abcd coeffs

[imagePoint]=findImagePoint(pointSource,n_ob,Hobj,Him,n_im,focalLength);

%check with estimated result
result2.imagePoint.position




%% PLOT the Entrance and Exit Pupil
EnP=result2.EntrancePupil;
ExP=result2.ExitPupil;

overDiam=[EnP.diam,ExP.diam];
maxDim=2*max(max(overDiam)); %twice the max pupil diameter

wave0=550; %nm   select a wavelengt

% set text box
textF='true';

% lens.draw
ppsfCamera.draw

% Graph pupil
[out1]=drawPupil(EnP,ExP,wave0,wave,textF);
%Point source
[out2]=drawPoint(pointSource,[],[],'y');
% Limiting rays
[out3]=drawComaRay(pointSource,EnP,wave0,wave,'y');
% and then the image point
[out4]=drawPoint(imagePoint,wave0,wave,'y');
% Limiting rays
[out5]=drawComaRay(imagePoint,ExP,wave0,wave,'y');


%% EXAMPLE: Image a new point
pSource0=[100 100 -100000000];
imagePoint0=findImagePoint(pSource0,n_ob,Hobj,Him,n_im,focalLength);
