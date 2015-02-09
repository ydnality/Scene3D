function obj = recordOnFilm(obj, film, nLines)
% Records the rays on the film surface
%
%  rays.recordOnFilm(obj, film, nLines)
%
% obj:    A ray
% film:   A sensor - can be either a plane or a spherical surface
% nLines: number of lines to draw for the illustration

%% Parameters 
if (ieNotDefined('nLines')), nLines = 0; end


%% Calculate intersection point of all the rays at sensor
intersectZ = repmat(film.position(3), [size(obj.origin, 1) 1]);
intersectT = (intersectZ - obj.origin(:, 3))./obj.direction(:, 3);

oldOrigin  = obj.origin;
obj.origin = obj.origin + obj.direction .* repmat(intersectT, [1 3]);


%% Plot phase space and for the future - ray-raced lines
if (nLines > 0)
    %parameters for plotting from lens to sensor
    lWidth = 0.1; lColor = [0 0.5 1]; lStyle = '-';

    obj.plotPhaseSpace();
    
    % These are the random sample of rays we will draw
    samps = obj.drawSamples;
    
    % These are the rays from the old origin to the new origin accounting
    % for the direction and current endpoint
    xCoordVector = [oldOrigin(samps,3) obj.origin(samps,3) NaN([nLines 1])]';
    yCoordVector = [oldOrigin(samps,2) obj.origin(samps,2) NaN([nLines 1])]';
    
    xCoordVector = real(xCoordVector(:));
    yCoordVector = real(yCoordVector(:));
    
    if isempty(obj.plotHandle), obj.plotHandle = vcNewGraphWin; end
    figure(obj.plotHandle);
    line(xCoordVector,  yCoordVector ,'Color',lColor,'LineWidth',lWidth,'LineStyle',lStyle);
    pause(0.1);
end

%% Remove dead rays
liveRays = obj.get('live rays');

% Record the real rays - if there are any
if(~isempty(liveRays.origin))
    if(isa(film, 'filmSphericalC'))
        
        % Project spherical sensor
        sensorRadius = film.radius;
        sensorCenter = film.get('sphericalCenter');  
        intersectPosition = liveRays.sphereIntersect(sensorCenter,sensorRadius);
        
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
        
        % When it is the plane, this ....
        intersectPosition = liveRays.origin;
        
        %imagePixel is the pixel that will gain a photon due to the traced ray
        imagePixel.position = [intersectPosition(:,2) intersectPosition(:, 1)];
        imagePixel.position = real(imagePixel.position); %add error handling for this
    else
        error('Invalid film type detected.  Quitting');
    end
    
    % Either way, we collect and count the photons incident on the film
    
    % 
    % This line takes a raw dimension and converts it to a position in
    % terms of pixels
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