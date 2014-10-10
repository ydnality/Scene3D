function [x_im,y_im] = psfEstimateImageCoord(coordType,P1,P2,wave,NA,exitPupilDiameter)
% Estimate image coordinates for the given pupil coordinate 
%
%  [x_im,y_im] = psfEstimateImageCoord(coordType,P1,P2,wave,NA,exitPupilDiameter)
%
% Not all the important cases are specified.
% The exit pupil case isn't properly built yet.
%
% INPUT
%   coordType: specify the type of input pupil coordinate.  The options are
%       'normalized cartesian'
%       'normalized polar'
%     P1: first coordinate of pupil function   {ro    or x}
%     P2: second coordinate of pupil function  {theta or y}
%     wave: wavelength vector
%     NA: Numerical apertures for different wavelengths ( column vector)
%     exitPupilDiameter: 
%        {1}: Exit Pupil Diameter (wavelength dependent) [needed for
%        convert real to normalized pupil coordinate] 
%
% OUTPUT
%  x_im: wavelength dependent x-coordinate in image space 
%  y_im: wavelength dependent y-coordinate in image space
%
% Examples:
%
%
% See also:  psfPupil2PSF, ComputeFFTpsf.m
%
% MP Vistasoft Team, 2014

%% Get some relevant parameter

NA = NA(:);
nW = size(NA,1);   % number of wavelengths

%% Main calculation

% Depending on the input coordinate type
coordType = ieParamFormat(coordType);

% The function names here don't make sense to me. Let's make it more clear.
% Mainly, psfNormalized2RealXXX - is it normalized in what size?  Real in
% what sense?  (BW).
switch coordType
    case{'normalizedcartesian','normcart'}
        % Does 'real' mean cartesian?? (BW)
        %Convert 'real' image coordinates to polar coordinates
        x_p = P1; y_p = P2; 
        % For each wavelength
         for li=1:nW       
            if (size(P1,1)==1) || (size(P1,2)==1) %case of INPUT as vector
                [x_im(li,:),y_im(li,:)] = psfNormalized2RealCoordinate(x_p,y_p,wave(li,1),NA(li,1));
            else % case of INPUT as matrix
                [x_im(:,:,li),y_im(:,:,li)] = psfNormalized2RealCoordinate(x_p,y_p,wave(li,1),NA(li,1));
            end
        end
    case{'normalizedpolar','normpol'}
        % Convert polar to cartesian coordinate
        ro=P1; theta=P2; 
        x_p=ro.*cos(theta);y_p=ro.*sin(theta);
        
        %Estimate 'real' image coordinate - Cartesian is real?
        for li=1:nW       
            if (size(P1,1)==1) || (size(P1,2)==1) %case of INPUT as vector
                [x_im(li,:),y_im(li,:)] = psfNormalized2RealCoordinate(x_p,y_p,wave(li,1),NA(li,1));
            else % case of INPUT as matrix
                [x_im(:,:,li),y_im(:,:,li)] = psfNormalized2RealCoordinate(x_p,y_p,wave(li,1),NA(li,1));
            end
        end

    case{'realcartesian','realcart'}
        % I think this one will require the exit pupil diameter.
        % And I don't understand these.
        error('Real Cartesian is not yet implemented');
    case{'realpolar'}
        % Or these.
        error('Real polar is not yet implemented');
    otherwise
        error ('Unknown coordinate type %s\n',coordType)
end


end



