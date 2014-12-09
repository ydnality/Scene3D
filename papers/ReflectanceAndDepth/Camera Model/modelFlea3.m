standardWave = 400:10:700;
standardWave = standardWave(:);

nWaves = length(standardWave);
deltaL = standardWave(2) - standardWave(1);

% Here are some of the key pixel properties of the PointGrey Flea3 camera
% with an ONSemi VITA 1300 1.3 Megapixel sensor.
wellCapacity   = 13700;                     % Electrons
conversiongain = 90*1e-6;                   % Volts/electron 
voltageSwing   = conversiongain*wellCapacity;
fillfactor     = 0.5;                       % This is a made up number
pixelSize      = 4.8*1e-6;                  % Meters
darkvoltage    = conversiongain*4.5;        % Volts/sec
readnoise      = 0.00096;                   % Volts
rows = 1024;                                 
cols = 1280; 
simRows = 105;                                 
simCols = 128;
exposureDuration = 0.06;                    % commented because we set autoexposure
dsnu =  conversiongain*30;                  % Volts (dark signal non-uniformity)
prnu = 2;                                   % Percent (ranging between 0 and 100) photodetector response non-uniformity
quantizationMethod = '8bit';
analogGain   = 1;                   % Used to adjust ISO speed
analogOffset = 0;                   % Used to account for sensor black level

noiseFlag = 2;                      % Use all the noise there is

qe = 0.53*ieReadSpectra(sprintf('%s/Camera Model/Responsivity',ReflDepthRootPath),standardWave);

% Some parameters of the lens used
focalLength = 0.07;
fNumber = 4;

% Multispectral properties
nChannels = 14;
channelCenterWaves = [447,470,505,530,590,627,655,680,780,850,880,940,395,365];

% Filter properties
nFilters = 8;
filterCenterWaves = {'',450,500,550,600,650,700,800};

cameraFilters = ones(nWaves,nFilters);
cameraFilters(:,1) = qe;
for f=2:nFilters
    fName = fullfile(ReflDepthRootPath,'Camera Model',sprintf('Bandpass%inm',filterCenterWaves{f}));
    cameraFilters(:,f) = diag(qe)*ieReadColorFilter(standardWave,fName);
end
