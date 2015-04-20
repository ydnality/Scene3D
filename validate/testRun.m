%% Apply proper PSF for every pixel

%dimensions of scene data
numRows = sceneGet(scene, 'rows');
numCols = sceneGet(scene,'cols');
%size of PSF film pixel
smallPixelSize = PSFStructure.film.size(1)/size(PSFStructure.film.image,1);     %mm
%size of scene sample pixel
largePixelSize = sceneGet(scene, 'sampleSize') * 10^3;  %mm
%amount that PSF film needs to be scaled to be equivalent to scene sample
%size
scaleFactor = smallPixelSize/largePixelSize;

%create an oi
oi = oiCreate;
oi = initDefaultSpectrum(oi);

%use scene data as photons first
photons = sceneGet(scene, 'photons');
oi = oiSet(oi, 'cphotons', photons);
padAmount = (round(size(PSFStructure.film.image,1) * scaleFactor) * 2);
oi = oiPad(oi, [padAmount padAmount]);
unBlurredPhotons = oiGet(oi, 'photons');
photonSum = zeros(size(oiGet(oi,'photons')));

vcAddObject(oi); oiWindow;

%loop through wavelength
for waveInd = 1:1   %length(wave)
    curWave = wave(waveInd);
    waveInd
    %loop through rows and cols
    for ii = 1:numRows
        for jj = 1:numCols

            %for each pixel, 
            %-figure out the field height(in angles)
            
            %-figure out the PSF for that field position and that depth
            %-reposition the PSF at that pixel and multiply it by
            %the pixel value.
            %-add to the sum
            
            curFieldHeightAngle = 0;   %fill this in in a bit
            
            depth = depthMap(ii, jj);
            
            %calculates the correct PSF for the current field angle, depth,
            %wavelength, given the PSf structure
            currentPSF = s3dPSFLookUp(curFieldHeightAngle, depth, curWave,PSFStructure);
            
            %scale PSF to the right size
            scaledCPSF = imresize(currentPSF, scaleFactor);
            scaledCPSF = scaledCPSF./sum(scaledCPSF(:)); %normalize
            scaledPSFNumRows = size(scaledCPSF, 1);
            
            
            %add to the sum
            centerPos = ceil(size(scaledCPSF, 1)/2);
            
            startingRow = ii + padAmount - centerPos;
            endingRow  = startingRow + scaledPSFNumRows-1;
            
            startingCol = jj + padAmount - centerPos;
            endingCol = startingCol + scaledPSFNumRows-1;
            
            photonSum(startingRow:endingRow,startingCol:endingCol, waveInd) = photonSum(startingRow:endingRow,startingCol:endingCol, waveInd)+ unBlurredPhotons(ii,jj, waveInd)*scaledCPSF;
        end
    end
end
    
oi = oiSet(oi,'cphotons',photonSum);
vcAddObject(oi); oiWindow;

%PSF lookup (can be separate function)

%figure out the PSF at that position and depth
%interpolate between depths and field positions by taking linear averages

%look up neighboring 4 PSF's (2x2)  2 field positions x 2 depths

%figure out where in that cube/square you are and interpolate

%return the interpolated value