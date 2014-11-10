function obj = recordOnFilm(obj, film, nLines)
%records the ray on the film surface
%nLines: number of lines to draw for the illustration

if (ieNotDefined('nLines'))
    nLines = 0;
end

%parameters for plotting from lens to sensor
lWidth = 0.1; lColor = [0 0.5 1]; lStyle = '-';



%make a clone of the rays
liveRays = rayC();
liveRays.makeDeepCopy(obj);

%calculate intersection point at sensor
intersectZ = repmat(film.position(3), [size(liveRays.origin, 1) 1]);
intersectT = (intersectZ - liveRays.origin(:, 3))./liveRays.direction(:, 3);
liveRays.origin = liveRays.origin + liveRays.direction .* repmat(intersectT, [1 3]);

%plot phase space and for the future - ray-raced lines
if (nLines > 0)
    liveRays.plotPhaseSpace();
    
    %draw debug lines
    % draw lines...TODO: figure this out...
    samps = obj.drawSamples;
    
    xCoordVector = [obj.origin(samps,3) liveRays.origin(samps,3) NaN([nLines 1])]';
    yCoordVector = [obj.origin(samps,2) liveRays.origin(samps,2) NaN([nLines 1])]';
    xCoordVector = real(xCoordVector(:));
    yCoordVector = real(yCoordVector(:));
    
    if isempty(obj.plotHandle), obj.plotHandle = vcNewGraphWin; end
    figure(obj.plotHandle);
    line(xCoordVector,  yCoordVector ,'Color',lColor,'LineWidth',lWidth,'LineStyle',lStyle);
    pause(0.1);
end

%remove dead rays
deadIndices = isnan(obj.waveIndex);
liveRays.origin(deadIndices, : ) = [];
liveRays.direction(deadIndices, : ) = [];
%liveRays.wavelength(deadIndices) = [];
liveRays.waveIndex(deadIndices) = [];

intersectPosition = liveRays.origin;

% Record the real rays - if there are any
if(~isempty(liveRays.origin))
    %imagePixel is the pixel that will gain a photon due to the traced ray
    imagePixel.position = [intersectPosition(:,2) intersectPosition(:, 1)];
    imagePixel.position = real(imagePixel.position); %add error handling for this
    
    %this line looks suspicious - what exactly does it do? will this cause
    %problems?
    imagePixel.position = round(imagePixel.position * film.resolution(2)/film.size(2) + ...
        repmat(-film.position(2:-1:1)*film.resolution(2)/film.size(2)  + (film.resolution(2:-1:1) + 1)./2, [size(imagePixel.position,1) 1]));   %
    
    imagePixel.wavelength = liveRays.get('wavelength');
    
    convertChannel = liveRays.waveIndex;
    %wantedPixel is the pixel that we wish to add 1 photon to
    wantedPixel = [imagePixel.position(:, 1) imagePixel.position(:,2) convertChannel];  %pixel to update
    
    recordablePixels =and(and(and(wantedPixel(:, 1) >= 1,  wantedPixel(:,1) <= film.resolution(1)), (wantedPixel(:, 2) > 1)), wantedPixel(:, 2) <= film.resolution(2));
    
    %remove the nonrecordablePixels
    wantedPixel = wantedPixel(recordablePixels, :);
    
    %correct for y coordinates
    wantedPixel(:, 1) =  film.resolution(1) + 1 - wantedPixel( :, 1);
    
    %make a histogram of wantedPixel in anticipation of adding
    %to film
    %  [count bins] = hist(single(wantedPixel));
    %  [count bins] = hist(single(wantedPixel), unique(single(wantedPixel), 'rows'));
    uniqueEntries =  unique(single(wantedPixel), 'rows');
    
    % Serializes the unique entries
    % I had an error here once where the wavelength in the
    % uniqueEntries was 7 but there were only 4 wavelength
    % dimensions in the film image.  Hmm. (BW).
    serialUniqueIndex = sub2ind(size(film.image), uniqueEntries(:,1), uniqueEntries(:,2), uniqueEntries(:,3));
    
    serialUniqueIndex = sort(serialUniqueIndex);
    
    serialWantedPixel = sub2ind(size(film.image), single(wantedPixel(:,1)), single(wantedPixel(:,2)), single(wantedPixel(:,3)));
    
    %special case for length 1.  For some reason, hist has
    %issues with length 1.
    if (length(serialUniqueIndex(:)) == 1)
        serializeFilm = film.image(:);
        %when there is only 1 bin, it doesn't matter how many photons, so just add one
        serializeFilm(serialUniqueIndex) = serializeFilm(serialUniqueIndex) + 1;
        film.image = reshape(serializeFilm, size(film.image));
    elseif(length(serialUniqueIndex(:) > 0))
        [countEntries] = hist(serialWantedPixel, serialUniqueIndex);
        %serialize the film, then the indices, then add by countEntries
        serializeFilm = film.image(:);
        serializeFilm(serialUniqueIndex) = serializeFilm(serialUniqueIndex) + countEntries';  %this line might be problematic...
        film.image = reshape(serializeFilm, size(film.image));
    else
        warning('No photons were collected on film!');
    end
end
end