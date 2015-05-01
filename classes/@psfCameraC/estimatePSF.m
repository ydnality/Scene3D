function estimatePSF(obj,nLines, jitterFlag, subsection, diffractionMethod)
% Estimate the PSF of a psfCamera
%
%   psfCamera.estimatePSF(obj)
%
% The camera has a point source, lens, and film.  This function saves the
% estimated psf into the psf camera object.  
%
% You can visualize it using the optical image derived from the psfCamera.
% To obtain the optical image of the PSF, use the following command after 
% calling this function:
%
%   psfCamera.oiCreate();
%
% AL/BW Vistasoft Team, Copyright 2014

if ieNotDefined('nLines'),     nLines = false;     end
if ieNotDefined('jitterFlag'), jitterFlag = false; end
if ieNotDefined('subsection'), subsection = []; end
if ieNotDefined('diffractionMethod'), diffractionMethod = 'HURB'; end

if (isequal(diffractionMethod, 'huygens') && obj.lens.diffractionEnabled) 
    % Trace from the point source to the entrance aperture of the
    % multielement lens
    
    lensMode = true;  %set false if ideal lens focused at infinity
    
    ppsfCFlag = false;
    obj.lens.diffractionEnabled = false;
    obj.rays = obj.lens.rtSourceToEntrance(obj.pointSource, ppsfCFlag, jitterFlag,[], subsection);

    % Duplicate the existing rays for each wavelength
    % Note that both lens and film have a wave, sigh.
    % obj.rays.expandWavelengths(obj.film.wave);
    obj.rays.expandWavelengths(obj.lens.wave);

    %lens intersection and raytrace
    obj.lens.rtThroughLens(obj.rays, nLines);

    % Something like this?  Extend the rays to the film plane?
    % if nLines > 0; obj.rays.draw(obj.film); end

    % intersect with "film" and add to film
    %obj.rays.recordOnFilm(obj.film, nLines);
    
    
    
    % Huygens ray-trace portion (PUT THIS IN A FUNCTION)
    % use a preset sensor size and pitch for now... somehow integrate this
    % with PSF camera later
    
    lambda = 550;  %for now
    
    %binSize   = [40000 40000];   %25;                 % nm
    binSize = [obj.film.size(1)/obj.film.resolution(1) obj.film.size(2)/obj.film.resolution(2)] .* 10^6;
    %numPixels = [50 50];                % In the sensor
    numPixels = [obj.film.resolution(1) obj.film.resolution(2)];
    numPixelsTot = numPixels(1) * numPixels(2);
    %imagePlaneDist     = 16*10^6; %100 * 10^6; %100mm
    imagePlaneDist = obj.film.position(3) * 10^6;
    %lensMode = false;


    
    %for each wavelength...
    
    for wIndex = 1:length(obj.lens.wave) 
        
        %estimated that the width of the 1st zero of airy disk will be .0336
        apXGridFlat = obj.rays.origin(:,1) * 10^6; %convert to nm
        apYGridFlat = obj.rays.origin(:,2) * 10^6;
        
        apXGridFlat = apXGridFlat(obj.rays.waveIndex == wIndex);
        apYGridFlat = apYGridFlat(obj.rays.waveIndex == wIndex);
        
        apXGridFlat = apXGridFlat(~isnan(apXGridFlat));
        apYGridFlat = apYGridFlat(~isnan(apYGridFlat));
        
        lambda = obj.lens.wave(wIndex);
        
        
        numApertureSamplesTot = length(apXGridFlat);

        %figure; plot(apXGridFlat, apYGridFlat, 'o');  %plot the aperture samples

        %create sensor grid
        % These are the locations on the sensor
        endLocations1DX = linspace(-numPixels(1)/2 * binSize(1), numPixels(1)/2 * binSize(1), numPixels(1));
        endLocations1DY = linspace(-numPixels(2)/2 * binSize(2), numPixels(2)/2 * binSize(2), numPixels(2));
        [endLGridX endLGridY] = meshgrid(endLocations1DX, endLocations1DY);
        endLGridXFlat = endLGridX(:);    %flatten the meshgrid
        endLGridYFlat = endLGridY(:);
        endLGridZFlat = ones(size(endLGridYFlat)) * imagePlaneDist;

        intensity = zeros(numPixels(1), numPixels(2), length(obj.lens.get('wave')));
        %intensityFlat = zeros(numPixels);
        %intensityFlat = intensityFlat(:);


        if (lensMode)
            initialD = obj.rays.distance(~isnan(obj.rays.distance));
        else
            initialD = zeros(numApertureSamplesTot, 1);
        end


        tic
        intensityFlat = zeros(numPixels);
        intensityFlat = intensityFlat(:);
        jobInterval = 5000;
        numJobs = ceil(numApertureSamplesTot/jobInterval);

        %split aperture into segments and do bsxfun on that, then combine later.
        %This can be parallelized later
        for job = 1:numJobs
            disp(job/numJobs)
            apXGridFlatCur = apXGridFlat((job-1) * jobInterval + 1:min(job * jobInterval, numApertureSamplesTot));
            apYGridFlatCur = apYGridFlat((job-1) * jobInterval + 1:min(job * jobInterval, numApertureSamplesTot));


            xDiffMat = bsxfun(@minus, endLGridXFlat, apXGridFlatCur');
            yDiffMat = bsxfun(@minus, endLGridYFlat, apYGridFlatCur');
            zDiffMat = repmat(endLGridZFlat, [1, length(apXGridFlatCur)]);

            initialDCur = initialD(((job-1) * jobInterval + 1:min(job * jobInterval, numApertureSamplesTot)));
            initialDMat = repmat(initialDCur', [length(endLGridXFlat), 1]);

            expD = exp(2 * pi * 1i .* ((sqrt(xDiffMat.^2 + yDiffMat.^2 + zDiffMat.^2) + initialDMat * 10^6)/lambda));
            intensityFlat = sum(expD, 2) + intensityFlat;
        end

        intensityFlat = abs(intensityFlat);
        toc

%         %plot results
        intensity1Wave = reshape(intensityFlat, size(intensity, 1), size(intensity, 2));
        obj.film.image(:,:,wIndex) = intensity1Wave;
% 
%         figure; imagesc(sqrt(intensity1Wave));
%         colormap('gray');
% 
%         if(lensMode)
%             lensModeText = 'Lens'
%         else
%             lensModeText = 'No Lens'
%         end
%         title([lensModeText ', ' int2str(imagePlaneDist/10^6) 'mm,' ' apDiameter: ' num2str(obj.lens.get('apertureMiddleD')) 'mm, pointSourceLocation: ' num2str(-obj.pointSource(3)) 'mm']);
%         xlabel([num2str(binSize(1) *numPixels(1)/10^6) 'mm']);
    end
    
    obj.lens.diffractionEnabled = true;
    %End Huygens ray-trace
else
    
    % Trace from the point source to the entrance aperture of the
    % multielement lens
    ppsfCFlag = false;
    % (pointSource, ppsfCFlag, jitterFlag, rtType, subSection, depthTriangles)
    obj.rays = obj.lens.rtSourceToEntrance(obj.pointSource, ppsfCFlag, jitterFlag, [], subsection);

    % Duplicate the existing rays for each wavelength
    % Note that both lens and film have a wave, sigh.
    % obj.rays.expandWavelengths(obj.film.wave);
    obj.rays.expandWavelengths(obj.lens.wave);

    %lens intersection and raytrace
    obj.lens.rtThroughLens(obj.rays, nLines);

    % Something like this?  Extend the rays to the film plane?
    % if nLines > 0; obj.rays.draw(obj.film); end

    % intersect with "film" and add to film
    obj.rays.recordOnFilm(obj.film, nLines);    
end
end