function PSF = s3dPSFLookUp(fieldHeight, depth, wavelength, PSFStructure)
%calculates a PSF given a PSF structure, field Height, depth, and
%wavelength
%
% What to do when we are outside the range?
%
% Also, this is the dopey interpolation.  Is there a better one that
% doesn't go all the way to VoLT
%
% 




% PSFStructure.fHAngle = PSFFieldHeightSamples;
% PSFStructure.depth = PSFDepthSamples;

%% figure out "anchor" depths

numDepths = length(PSFStructure.depth);
for i = 1:(numDepths-1)
    d = PSFStructure.depth(i);
    nextD = PSFStructure.depth(i + 1);
    if (d <= depth && depth < nextD)
       break; 
    end
end

if (depth <= PSFStructure.depth(1))
    lesserDepthInd = 1;
    greaterDepthInd = 1;
elseif (depth >= PSFStructure.depth(numDepths))
     lesserDepthInd = numDepths;
     greaterDepthInd = numDepths;
else
    lesserDepthInd = i;
    greaterDepthInd = i+1;
end


%% figure out "anchor" field heigths
numHeights= length(PSFStructure.fHAngle);

for i = 1:(numHeights-1)
    fH = PSFStructure.fHAngle(i);
    nextFH = PSFStructure.fHAngle(i+1);
    if (fH <= fieldHeight && fieldHeight < nextFH)
       break; 
    end
end

if (fieldHeight <= PSFStructure.fHAngle(1))
    lesserHeightInd = 1;
    greaterHeightInd = 1;
elseif (fieldHeight >= PSFStructure.fHAngle(numHeights))
     lesserHeightInd = numHeights;
     greaterHeightInd = numHeights;
else
    lesserHeightInd = i;
    greaterHeightInd = i+1;
end

%% figure out "anchor" wavelength - TODO: make this into function later
numWave= length(PSFStructure.wave);

for i = 1:(numWave-1)
    wave = PSFStructure.wave(i);
    nextWave = PSFStructure.wave(i+1);
    if (wave <= wavelength && wavelength < nextWave)
       break; 
    end
end

if (wavelength <= PSFStructure.wave(1))
    lesserWaveInd = 1;
    greaterWaveInd = 1;
elseif (wavelength >= PSFStructure.wave(numWave))
     lesserWaveInd = numWave;
     greaterWaveInd = numWave;
else
    lesserWaveInd = i;
    greaterWaveInd = i+1;
end


%% linear interpolation
depthStep =  PSFStructure.depth(greaterDepthInd) - PSFStructure.depth(lesserDepthInd);
if (depthStep ==0)
    lesserDepthWeight = 1;
    greaterDepthWeight = 0;
else
    lesserDepthWeight = (PSFStructure.depth(greaterDepthInd) - depth)/depthStep;  
    greaterDepthWeight = (depth - PSFStructure.depth(lesserDepthInd))/depthStep;
end

heightStep =  PSFStructure.fHAngle(greaterHeightInd) - PSFStructure.fHAngle(lesserHeightInd);
if(heightStep ==0) 
    lesserHeightWeight = 1;
    greaterHeightWeight = 0;
else
    lesserHeightWeight = (PSFStructure.fHAngle(greaterHeightInd) - fieldHeight)/heightStep;  
    greaterHeightWeight = (fieldHeight - PSFStructure.fHAngle(lesserHeightInd))/heightStep;
end

waveStep = PSFStructure.wave(greaterWaveInd) - PSFStructure.wave(lesserWaveInd);
if(waveStep ==0) 
    lesserWaveWeight = 1;
    greaterWaveWeight = 0;
else
    lesserWaveWeight = (PSFStructure.wave(greaterWaveInd) - wavelength)/waveStep;  
    greaterWaveWeight = (wavelength - PSFStructure.wave(lesserWaveInd))/waveStep;
end


%  PSF(:,:,waveInd,depthInd, fHIndex) = curPSF;

%ignore wavelength for now.  TODO: Add it in later
PSF = lesserWaveWeight *    (PSFStructure.PSF( :,:, lesserWaveInd ,lesserDepthInd, lesserHeightInd) * lesserDepthWeight * lesserHeightWeight + ...
                              PSFStructure.PSF( :,:, lesserWaveInd ,lesserDepthInd, greaterHeightInd) * lesserDepthWeight * greaterHeightWeight + ...
                              PSFStructure.PSF( :,:, lesserWaveInd ,greaterDepthInd, lesserHeightInd) * greaterDepthWeight * lesserHeightWeight + ...
                              PSFStructure.PSF( :,:, lesserWaveInd ,greaterDepthInd, greaterHeightInd) * greaterDepthWeight * greaterHeightWeight) + ...
      greaterWaveWeight * (PSFStructure.PSF( :,:, greaterWaveInd ,lesserDepthInd, lesserHeightInd) * lesserDepthWeight * lesserHeightWeight + ...
                              PSFStructure.PSF( :,:, greaterWaveInd ,lesserDepthInd, greaterHeightInd) * lesserDepthWeight * greaterHeightWeight + ...
                              PSFStructure.PSF( :,:, greaterWaveInd ,greaterDepthInd, lesserHeightInd) * greaterDepthWeight * lesserHeightWeight + ...
                              PSFStructure.PSF( :,:, greaterWaveInd ,greaterDepthInd, greaterHeightInd) * greaterDepthWeight * greaterHeightWeight);
  
%normalize to sum to 1
PSF = PSF./sum(PSF(:));
 
end