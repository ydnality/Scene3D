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
    if(isa(film, 'filmSphericalC'))
        
        %%
        % Project sphereSense

        % sensor properties
        
        sensorRadius = film.radius;
        sensorCenter = film.get('sphericalCenter');  %TODO: have this be a get?
        
        % intersect with a spherical surface and find intersection point
        
        % Spherical element
        %repCenter = repmat(curEl.sCenter, [nRays 1]);
        nRays = length(liveRays.origin);
        repCenter = repmat(sensorCenter, [nRays 1]);
        
        repRadius = repmat(sensorRadius, [nRays 1]);
        
        % Radicand from vector form of Snell's Law
        radicand = dot(liveRays.direction, liveRays.origin - repCenter, 2).^2 - ...
            ( dot(liveRays.origin - repCenter, liveRays.origin -repCenter, 2)) + repRadius.^2;
        
        % Calculate something about the ray angle with respect
        % to the current surface.  AL to figure this one out
        % and put in a book reference.
        if (sensorRadius < 0)
            intersectT = (-dot(liveRays.direction, liveRays.origin - repCenter, 2) + sqrt(radicand));
        else
            intersectT = (-dot(liveRays.direction, liveRays.origin - repCenter, 2) - sqrt(radicand));
        end
        
        %make sure that intersectT is > 0   This does not apply for this
        %case, because sometimes a curved sensor is in front of the flat
        %sensor plane
%         if (min(intersectT(:)) < 0)
%             fprintf('intersectT less than 0 for lens %i');
%         end
        
        repIntersectT = repmat(intersectT, [1 3]);
        intersectPosition = liveRays.origin + repIntersectT .* liveRays.direction;
        
        if(nLines > 0)
            figure(obj.plotHandle);    hold on; plot(intersectPosition(:,3), intersectPosition(:,2), '.');
        end
        % Convert intersection point into spherical coordinate system...
        x = intersectPosition(:,1);
        y = intersectPosition(:,2);
        z = intersectPosition(:,3) - sensorCenter(3);
        r = sqrt(x.^2 + y.^2 + z.^2);
        theta = atan(x./z);
        phi = asin(y./r);
        
        imagePixel.position = [phi * abs(sensorRadius) theta * abs(sensorRadius) ];
        figure; hist(theta);
        figure; hist(phi);
        
        % Record on sensor
        
    elseif(isa(film, 'filmC'))
        
        %imagePixel is the pixel that will gain a photon due to the traced ray
        imagePixel.position = [intersectPosition(:,2) intersectPosition(:, 1)];
        imagePixel.position = real(imagePixel.position); %add error handling for this
    else
        error('Invalid film type detected.  Quitting');
    end
    
    
    %this line takes a raw dimension and converts it to a position in
    %terms of pixels
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